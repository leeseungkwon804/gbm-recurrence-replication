#!/usr/bin/env python3
"""How far is tumour-adjacent normal from healthy tissue, relative to how far the
tumour is?

One pipeline for all three classes (Xena toil recompute), so quantification does
not differ. The honest caveat up front: GTEx is post-mortem and TCGA is surgical,
so any gap between adjacent normal and GTEx mixes field effect with collection
mode. This first pass only measures the magnitude -- if adjacent normal sits
essentially on top of healthy, the question dies here and cheaply.
"""
import pandas as pd, numpy as np, gzip, os
from scipy import stats

D = '/workspace/_scout_ksh/tn'
TUM, NAT, HEA = 'Primary Tumor', 'Solid Tissue Normal', 'Normal Tissue'
TISSUES = ['Breast', 'Lung', 'Prostate', 'Liver', 'Colon', 'Stomach', 'Kidney']

p = pd.read_csv(f'{D}/pheno.txt.gz', sep='\t', encoding='latin-1', low_memory=False)
p.columns = [c.strip('_') for c in p.columns]
p = p.rename(columns={'primary_site': 'site', 'sample_type': 'stype', 'study': 'study'})
sel = p[p.stype.isin([TUM, NAT, HEA]) & p.site.isin(TISSUES)]
want = set(sel['sample'].astype(str))
print(f'필요 샘플 {len(want):,}개 · 장기 {len(TISSUES)}개')

f = f'{D}/TcgaTargetGtex_rsem_gene_tpm.gz'
# header first, to pick columns
with gzip.open(f, 'rt') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
cols = [hdr[0]] + [c for c in hdr[1:] if c in want]
print(f'행렬에서 찾은 샘플 {len(cols)-1:,}개')

E = pd.read_csv(f, sep='\t', usecols=cols, index_col=0)
print(f'발현 {E.shape[0]:,} × {E.shape[1]:,}   (log2(TPM+0.001))')

meta = sel.set_index('sample').loc[[c for c in E.columns]]
# drop low-expression genes globally
keep = (E > 0).sum(axis=1) >= E.shape[1] * 0.5
E = E[keep]
print(f'발현 검출 유전자 {E.shape[0]:,}')
print()

print(f'{"장기":<10}{"종양":>6}{"인접정상":>8}{"건강":>6}'
      f'{"D(인접,건강)":>13}{"D(종양,건강)":>13}{"비율":>8}{"인접≠건강 유전자":>16}')
out = []
for t in TISSUES:
    m = meta[meta.site == t]
    a = E[m[m.stype == TUM].index]
    b = E[m[m.stype == NAT].index]
    c = E[m[m.stype == HEA].index]
    if b.shape[1] < 10 or c.shape[1] < 10:
        continue
    ca, cb, cc = a.mean(axis=1), b.mean(axis=1), c.mean(axis=1)
    d_nat = np.linalg.norm(cb - cc) / np.sqrt(len(cc))
    d_tum = np.linalg.norm(ca - cc) / np.sqrt(len(cc))
    # per-gene test: adjacent normal vs healthy
    tt, pp = stats.ttest_ind(b.values, c.values, axis=1, equal_var=False)
    from statsmodels.stats.multitest import multipletests
    q = multipletests(np.nan_to_num(pp, nan=1.0), method='fdr_bh')[1]
    nsig = int((q < 0.01).sum())
    out.append((t, a.shape[1], b.shape[1], c.shape[1], d_nat, d_tum, d_nat / d_tum, nsig))
    print(f'{t:<10}{a.shape[1]:>6}{b.shape[1]:>8}{c.shape[1]:>6}'
          f'{d_nat:>13.3f}{d_tum:>13.3f}{d_nat/d_tum:>8.2f}{nsig:>16,}')

r = pd.DataFrame(out, columns=['tissue', 'n_tum', 'n_nat', 'n_hea',
                               'd_nat_hea', 'd_tum_hea', 'ratio', 'n_sig_q01'])
r.to_csv(f'{D}/distance_summary.tsv', sep='\t', index=False)
print()
print(f'  → 인접정상이 건강에서 떨어진 정도가 종양의 {r.ratio.min():.0%}~{r.ratio.max():.0%}')
print(f'  → 인접정상 대 건강에서 q<0.01 유전자 중앙값 {r.n_sig_q01.median():,.0f} '
      f'(검정 유전자 {E.shape[0]:,})')
