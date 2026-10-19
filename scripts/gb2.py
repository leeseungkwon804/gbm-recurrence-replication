#!/usr/bin/env python3
"""Is the first-to-second-resection change predictable in a new patient?

Learn the average within-patient change in one cohort, then ask, for each patient
in the other cohort -- different people, different platform -- how well that
learned change matches their own. Nothing from the test cohort enters the
training, so there is no leakage, and correlation is used so the two platforms'
different scales do not matter.

A gene-label permutation gives the null.
"""
import numpy as np, pandas as pd, json
from scipy import stats

D = '/workspace/_scout_ksh'
dA = pd.read_csv(f'{D}/out/A_deltas.tsv', sep='\t', index_col=0)   # ensembl x patients
dB = pd.read_csv(f'{D}/out/B_deltas.tsv', sep='\t', index_col=0)   # symbol x patients
sym = json.load(open(f'{D}/out/ens2sym.json'))
dA.index = [sym.get(g.split('.')[0], '') for g in dA.index]
dA = dA[dA.index != ''].groupby(level=0).mean()
shared = sorted(set(dA.index) & set(dB.index))
A, B = dA.loc[shared], dB.loc[shared]
print(f'공통 유전자 {len(shared)} · A {A.shape[1]}쌍 · B {B.shape[1]}쌍')
print()


def crosstest(train, test, ltr, lte, rng):
    d = train.mean(axis=1).values                      # learned change, training cohort only
    rows = []
    for p in test.columns:
        v = test[p].values
        r, _ = stats.spearmanr(d, v)
        rows.append(r)
    rows = np.array(rows)
    # permutation null: shuffle the learned vector's gene labels
    null = []
    for _ in range(2000):
        dp = rng.permutation(d)
        null.append(np.mean([stats.spearmanr(dp, test[p].values)[0] for p in test.columns[:4]]))
    null = np.array(null)
    w = stats.wilcoxon(rows, alternative='greater')
    print(f'=== {ltr} 에서 배운 변화로 {lte} 환자를 맞추나 ===')
    print(f'    환자별 상관 중앙 {np.median(rows):+.3f} · 범위 {rows.min():+.3f}~{rows.max():+.3f}')
    print(f'    양의 상관 환자 {int((rows>0).sum())}/{len(rows)}')
    print(f'    부호순위검정 p = {w.pvalue:.3g}')
    print(f'    순열 귀무 중앙 {np.median(null):+.3f} (95% {np.percentile(null,2.5):+.3f}~'
          f'{np.percentile(null,97.5):+.3f})')
    return rows


rng = np.random.default_rng(0)
rAB = crosstest(A, B, 'A(RNA-seq 14쌍)', 'B(NanoString 17쌍)', rng)
print()
rBA = crosstest(B, A, 'B(NanoString 17쌍)', 'A(RNA-seq 14쌍)', rng)
print()

# 얼마나 많은 환자가 "평균 변화"로 설명되나
print('=== 개인차는 얼마나 되나 ===')
for lab, M in [('A', A), ('B', B)]:
    m = M.mean(axis=1).values
    rs = [stats.spearmanr(m, M[p].values)[0] for p in M.columns]
    print(f'  {lab}: 자기 코호트 평균과의 상관 중앙 {np.median(rs):.3f} '
          f'(범위 {min(rs):.3f}~{max(rs):.3f})')
print()
print('  → 교차코호트 상관이 자기코호트 상관에 가까울수록,')
print('    변화가 코호트 특성이 아니라 질환 공통이라는 뜻이다')
