#!/usr/bin/env python3
"""Which tissues support the three-way comparison, and is the TCGA tumour/normal
pairing recoverable?

Three sample classes, one pipeline: TCGA tumour, TCGA adjacent normal, GTEx
healthy. A tissue is usable only if adjacent normal and healthy are both present
in workable numbers. TCGA barcodes carry the patient in the first twelve
characters, which is what makes the tumour-versus-its-own-adjacent-normal
contrast possible.
"""
import pandas as pd, numpy as np, re

p = pd.read_csv('/workspace/_scout_ksh/tn/pheno.txt.gz', sep='\t',
                encoding='latin-1', low_memory=False)
p.columns = [c.strip('_') for c in p.columns]
p = p.rename(columns={'primary disease or tissue': 'disease', 'primary_site': 'site',
                      'sample_type': 'stype', 'study': 'study'})

TUM = 'Primary Tumor'
NAT = 'Solid Tissue Normal'
HEA = 'Normal Tissue'
sel = p[p.stype.isin([TUM, NAT, HEA])].copy()

tab = sel.pivot_table(index='site', columns='stype', values='sample',
                      aggfunc='count').fillna(0).astype(int)
tab = tab[[TUM, NAT, HEA]]
tab['min3'] = tab.min(axis=1)
tab = tab.sort_values('min3', ascending=False)
print('=== 장기별 3자 표본수 (종양 / TCGA 인접정상 / GTEx 건강) ===')
print(tab[tab.min3 >= 5].to_string())
print()
ok = tab[(tab[NAT] >= 15) & (tab[HEA] >= 20)]
print(f'3자 비교 가능 장기 (인접정상>=15, 건강>=20): {len(ok)}개')
print(f'  {list(ok.index)}')
print()

# TCGA 환자 짝지음
t = sel[(sel.study == 'TCGA') & sel.stype.isin([TUM, NAT])].copy()
t['patient'] = t['sample'].astype(str).str.slice(0, 12)
g = t.groupby(['site', 'patient']).stype.nunique()
paired = g[g == 2].reset_index().groupby('site').size().sort_values(ascending=False)
print('=== 같은 환자에서 종양+인접정상이 모두 있는 수 ===')
print(paired[paired >= 5].to_string())
print(f'\n  총 짝지은 환자 {int(paired.sum())}명 · 장기 {len(paired)}개')
