#!/usr/bin/env python3
"""Is the angiogenic axis systematically the least reproducible part of the
primary-to-recurrence change, or were VEGFA and ANGPT1 just two unlucky genes?

Two independent cohorts on two platforms, 699 genes in common. For each
pre-specified programme, measure how well the within-patient change agrees
between cohorts, then ask whether angiogenesis sits outside the rest.

Gene sets are canonical and fixed before looking at the per-set result.
"""
import numpy as np, pandas as pd
from scipy import stats

D = '/workspace/_scout_ksh'
j = pd.read_csv(f'{D}/out/replication.tsv', sep='\t', index_col=0)   # A_delta,A_q,B_delta,B_q

SETS = {
 'angiogenesis': ['VEGFA', 'VEGFB', 'VEGFC', 'ANGPT1', 'ANGPT2', 'KDR', 'FLT1', 'TEK',
                  'PECAM1', 'ESM1', 'CD34', 'ENG', 'NOTCH4', 'DLL4', 'HIF1A', 'EPAS1'],
 'myeloid': ['CSF1R', 'SIRPA', 'MS4A6A', 'IRF8', 'DAB2', 'CMKLR1', 'AIF1', 'ITGAM',
             'CD163', 'MRC1', 'TYROBP', 'FCGR3A', 'C1QA', 'C1QB', 'CD68', 'TREM2'],
 'T cell': ['CD3D', 'CD3E', 'CD8A', 'CD2', 'LCK', 'ZAP70', 'IL7R', 'CCL5', 'GZMA',
            'GZMB', 'PRF1', 'CTLA4', 'PDCD1', 'LAG3', 'TIGIT', 'FOXP3'],
 'stemness': ['SOX2', 'NOTCH1', 'OLIG2', 'NES', 'POU3F2', 'DLL3', 'ASCL1', 'PROM1',
              'MYC', 'BMI1', 'FUT4', 'L1CAM'],
 'ECM/mesenchymal': ['VCAN', 'LOXL2', 'CDH11', 'FN1', 'COL1A1', 'COL1A2', 'COL3A1',
                     'SERPINE1', 'TGFB1', 'CD44', 'TIMP1', 'MMP2', 'MMP9', 'SPARC'],
 'proliferation': ['MKI67', 'TOP2A', 'CCNB1', 'CDK1', 'AURKA', 'BUB1', 'PLK1', 'TYMS',
                   'RRM2', 'PCNA', 'TFRC', 'ANLN'],
 'interferon': ['STAT1', 'STAT2', 'IRF1', 'IRF7', 'ISG15', 'MX1', 'OAS1', 'IFIT1',
                'IFIT3', 'CXCL9', 'CXCL10', 'GBP1'],
 'hypoxia': ['CA9', 'SLC2A1', 'LDHA', 'PGK1', 'VEGFA', 'ADM', 'NDRG1', 'PDK1', 'BNIP3'],
}

rows = []
for nm, gs in SETS.items():
    g = [x for x in gs if x in j.index]
    if len(g) < 5:
        rows.append((nm, len(g), np.nan, np.nan, np.nan)); continue
    s = j.loc[g]
    agree = ((s.A_delta > 0) == (s.B_delta > 0)).mean()
    rho, _ = stats.spearmanr(s.A_delta, s.B_delta)
    # magnitude of disagreement, scaled
    disc = np.median(np.abs(s.A_delta - s.B_delta))
    rows.append((nm, len(g), agree, rho, disc))

r = pd.DataFrame(rows, columns=['set', 'n', 'agree', 'rho', 'median_gap']).set_index('set')
r = r.sort_values('agree')
print('=== 프로그램별 두 코호트 일치 (699개 공통 유전자 안에서) ===')
print(r.to_string(float_format=lambda v: f'{v:.3f}'))
print()

allagree = ((j.A_delta > 0) == (j.B_delta > 0)).mean()
print(f'전체 699개 방향 일치율 {allagree:.3f}')
print()
print('=== 혈관신생이 나머지와 다른가 (이항검정) ===')
ang = [x for x in SETS['angiogenesis'] if x in j.index]
oth = [x for x in j.index if x not in ang]
a_ok = int(((j.loc[ang].A_delta > 0) == (j.loc[ang].B_delta > 0)).sum())
o_ok = int(((j.loc[oth].A_delta > 0) == (j.loc[oth].B_delta > 0)).sum())
tab = [[a_ok, len(ang) - a_ok], [o_ok, len(oth) - o_ok]]
odds, p = stats.fisher_exact(tab)
print(f'  혈관신생 {a_ok}/{len(ang)} 일치 · 나머지 {o_ok}/{len(oth)} 일치')
print(f'  Fisher p = {p:.4g}  오즈비 {odds:.3f}')
print()
print('=== 혈관신생 유전자 하나씩 ===')
print(j.loc[ang].assign(같은방향=lambda x: np.where((x.A_delta > 0) == (x.B_delta > 0), 'O', 'X'))
      .to_string(float_format=lambda v: f'{v:.3f}'))
