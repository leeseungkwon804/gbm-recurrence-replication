#!/usr/bin/env python3
"""GTEx brain: check the power before the question.

Same donor, many brain regions, one collection event. That removes the batch
problem that killed the TCGA-versus-GTEx comparison earlier today, because the
regions being compared come out of the same body at the same time.

What matters first: how many donors actually contribute several regions, and what
donor-level variables are attached.
"""
import pandas as pd, numpy as np, collections

D = '/workspace/_scout_ksh/tn'
p = pd.read_csv(f'{D}/pheno.txt.gz', sep='\t', encoding='latin-1', low_memory=False)
p.columns = [c.strip('_') for c in p.columns]
p = p.rename(columns={'primary_site': 'site', 'sample_type': 'stype',
                      'primary disease or tissue': 'detail'})
g = p[(p.study == 'GTEX') & (p.site == 'Brain')].copy()
print(f'GTEx 뇌 샘플 {len(g):,}')
print('\n=== 영역별 샘플수 ===')
print(g.detail.value_counts().to_string())

# GTEx 샘플 ID: GTEX-<donor>-<...>
g['donor'] = g['sample'].astype(str).str.split('-').str[:2].str.join('-')
per = g.groupby('donor').detail.nunique().sort_values(ascending=False)
print(f'\n=== 공여자 {per.size}명 · 1인당 영역수 ===')
print(f'  중앙 {per.median():.0f} · 범위 {per.min()}-{per.max()}')
print(f'  영역 5개 이상 보유 공여자 {int((per >= 5).sum())}명')
print(f'  영역 10개 이상 보유 공여자 {int((per >= 10).sum())}명')

# 영역쌍별로 같은 공여자가 몇 명이나 되나
regions = list(g.detail.value_counts().index)
pair = collections.Counter()
byd = g.groupby('donor').detail.apply(set)
for s in byd:
    for i, a in enumerate(regions):
        for b in regions[i + 1:]:
            if a in s and b in s:
                pair[(a, b)] += 1
top = pair.most_common(10)
print('\n=== 같은 공여자에서 둘 다 있는 영역쌍 (상위 10) ===')
for (a, b), n in top:
    print(f'  {n:>4}명  {a[:34]:36s} × {b[:34]}')
print(f'\n  공여자 20명 이상 공유하는 영역쌍 {sum(1 for v in pair.values() if v >= 20)}개')
