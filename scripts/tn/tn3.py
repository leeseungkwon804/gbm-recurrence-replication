#!/usr/bin/env python3
"""Batch or biology?

A handling artefact pushes the same genes the same way in every organ, because
it is a property of how the sample was collected, not of the tissue. A field
effect around a tumour should look different in breast than in liver. So:
correlate the adjacent-normal-minus-healthy shift between organ pairs. High
correlation across unrelated organs means the signal is technical.

As a second read, compare the shift against the tumour-minus-healthy shift; if
adjacent normal is simply carrying the same artefact, the two will also track.
"""
import pandas as pd, numpy as np, gzip
from scipy import stats

D = '/workspace/_scout_ksh/tn'
TUM, NAT, HEA = 'Primary Tumor', 'Solid Tissue Normal', 'Normal Tissue'
TISSUES = ['Breast', 'Lung', 'Prostate', 'Liver', 'Colon', 'Stomach', 'Kidney']

p = pd.read_csv(f'{D}/pheno.txt.gz', sep='\t', encoding='latin-1', low_memory=False)
p.columns = [c.strip('_') for c in p.columns]
p = p.rename(columns={'primary_site': 'site', 'sample_type': 'stype'})
sel = p[p.stype.isin([TUM, NAT, HEA]) & p.site.isin(TISSUES)]
want = set(sel['sample'].astype(str))
f = f'{D}/TcgaTargetGtex_rsem_gene_tpm.gz'
with gzip.open(f, 'rt') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
E = pd.read_csv(f, sep='\t', usecols=[hdr[0]] + [c for c in hdr[1:] if c in want],
                index_col=0)
meta = sel.set_index('sample').loc[list(E.columns)]
E = E[(E > 0).sum(axis=1) >= E.shape[1] * 0.5]

shift_nat, shift_tum = {}, {}
for t in TISSUES:
    m = meta[meta.site == t]
    a = E[m[m.stype == TUM].index].mean(axis=1)
    b = E[m[m.stype == NAT].index].mean(axis=1)
    c = E[m[m.stype == HEA].index].mean(axis=1)
    shift_nat[t] = (b - c)
    shift_tum[t] = (a - c)
S = pd.DataFrame(shift_nat)
T = pd.DataFrame(shift_tum)

print('=== 인접정상−건강 변화가 장기 사이에서 얼마나 닮았나 (Spearman) ===')
print('    (배치면 높고, 장기별 생물학이면 낮아야 한다)')
C = S.corr(method='spearman')
print(C.round(2).to_string())
iu = np.triu_indices_from(C.values, k=1)
print(f'\n  장기쌍 상관 중앙값 {np.median(C.values[iu]):.3f} '
      f'(범위 {C.values[iu].min():.2f}~{C.values[iu].max():.2f})')
print()

print('=== 같은 장기 안에서, 인접정상 변화와 종양 변화가 닮았나 ===')
for t in TISSUES:
    r, _ = stats.spearmanr(S[t], T[t])
    print(f'  {t:<10} rho={r:.3f}')
print()

# 무관한 장기쌍 하나를 직접 본다
r_bl, _ = stats.spearmanr(S['Breast'], S['Liver'])
r_pk, _ = stats.spearmanr(S['Prostate'], S['Kidney'])
print(f'  유방 대 간 rho={r_bl:.3f} · 전립선 대 신장 rho={r_pk:.3f}')
print('  → 서로 무관한 장기인데도 높으면 조직 생물학이 아니라 채취·처리 신호다')
