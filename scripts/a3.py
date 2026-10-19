#!/usr/bin/env python3
"""Third pass: how reproducible is each change, and does it hold at module level."""
import json
import numpy as np, pandas as pd
from scipy import stats

D = '/workspace/_scout_ksh'
j = pd.read_csv(f'{D}/out/replication.tsv', sep='\t', index_col=0)

sig = j[j.B_q < 0.05]
same = (sig.A_delta > 0) == (sig.B_delta > 0)
k, n = int(same.sum()), len(sig)
print(f'=== B 에서 q<0.05 인 {n}개가 A 에서 같은 방향인가 ===')
print(f'   {k}/{n} = {100 * k / n:.1f}%  이항검정 p = {stats.binomtest(k, n, 0.5).pvalue:.3g}')
top = j.sort_values('B_q').head(30)
k2 = int(((top.A_delta > 0) == (top.B_delta > 0)).sum())
print(f'   상위 30: {k2}/30  p = {stats.binomtest(k2, 30, 0.5).pvalue:.3g}')
print()
print('=== 방향이 어긋난 것 ===')
print(sig[~same][['A_delta', 'A_q', 'B_delta', 'B_q']]
      .to_string(float_format=lambda v: f'{v:.3f}'))
print()

MYE = ['CSF1R', 'SIRPA', 'MS4A6A', 'IRF8', 'DAB2', 'CMKLR1', 'AIF1', 'ITGAM',
       'CD163', 'MRC1', 'TYROBP', 'FCGR3A', 'C1QA', 'C1QB']
STEM = ['SOX2', 'NOTCH1', 'OLIG2', 'NES', 'POU3F2', 'DLL3', 'ASCL1']
MES = ['VCAN', 'LOXL2', 'CDH11', 'FN1', 'COL1A1', 'SERPINE1', 'TGFB1', 'CD44']

dA = pd.read_csv(f'{D}/out/A_deltas.tsv', sep='\t', index_col=0)
m = json.load(open(f'{D}/out/ens2sym.json'))
dA['symbol'] = [m.get(g.split('.')[0], '') for g in dA.index]
dA = dA[dA.symbol != ''].groupby('symbol').mean()
dB = pd.read_csv(f'{D}/out/B_deltas.tsv', sep='\t', index_col=0)

print('=== 모듈 수준 (환자별 평균 Δ, 부호순위검정) ===')
for nm, gs in [('골수계/대식세포', MYE), ('줄기성', STEM), ('간엽/기질', MES)]:
    print(f'  {nm}')
    for lab, d in [('A RNA-seq  14쌍', dA), ('B NanoString 17쌍', dB)]:
        g = [x for x in gs if x in d.index]
        if len(g) < 3:
            print(f'     {lab}: 측정된 유전자 {len(g)}개뿐 — 건너뜀')
            continue
        v = d.loc[g].mean().dropna().values.astype(float)
        p = stats.wilcoxon(v, alternative='two-sided').pvalue
        print(f'     {lab}: Δ중앙 {np.median(v):+.3f} · 상승환자 {100 * (v > 0).mean():.0f}% '
              f'· p={p:.4f}  (유전자 {len(g)}개)')
