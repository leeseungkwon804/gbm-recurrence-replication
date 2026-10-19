#!/usr/bin/env python3
"""Second pass: replication across two platforms, on a pre-specified gene set.

Cohort A (GSE339484, RNA-seq, 14 pairs) cannot reach q < 0.05 on 23,084 genes:
the smallest p an exact signed-rank test can return at n = 14 is 1.22e-4, so the
smallest attainable q is 0.123. That is a power floor, not an absence of signal --
939 genes reach p < 0.01 where 231 are expected.

So the design is the one C2 used: fix the gene set before testing. Cohort B is
the NanoString PanCancer IO 360 panel, 770 genes chosen to read the immune
microenvironment. Restricting cohort A to those same genes drops the testing
burden about 29-fold and makes the two cohorts directly comparable on a set that
neither of us chose after seeing the data.

Cohort B also carries relapse time in days, which separates 'the tumour changed'
from 'time passed'.
"""
import os, json, urllib.request, urllib.parse
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

D = '/workspace/_scout_ksh'
OUT = f'{D}/out'
dA = pd.read_csv(f'{OUT}/A_deltas.tsv', sep='\t', index_col=0)
dB = pd.read_csv(f'{OUT}/B_deltas.tsv', sep='\t', index_col=0)
rB = pd.read_csv(f'{OUT}/B_GSE254874_pergene.tsv', sep='\t', index_col=0)


def symbols(ens):
    out = {}
    ens = [e.split('.')[0] for e in ens]
    for i in range(0, len(ens), 900):
        body = urllib.parse.urlencode({'ids': ','.join(ens[i:i + 900]),
                                       'fields': 'symbol', 'species': 'human'}).encode()
        try:
            for h in json.load(urllib.request.urlopen('https://mygene.info/v3/gene', body, timeout=120)):
                if h.get('symbol'):
                    out[h['query']] = h['symbol']
        except Exception as e:
            print('   mygene 실패:', e)
            break
    return out


cache = f'{OUT}/ens2sym.json'
if os.path.exists(cache):
    m = json.load(open(cache))
else:
    m = symbols(list(dA.index))
    json.dump(m, open(cache, 'w'))
print(f'Ensembl→symbol 매핑 {len(m):,}/{len(dA):,}')

dA = dA.copy()
dA['symbol'] = [m.get(g.split('.')[0], '') for g in dA.index]
dA = dA[dA.symbol != ''].groupby('symbol').mean()      # collapse duplicate symbols
panel = sorted(set(dB.index) & set(dA.index))
print(f'패널 {len(dB)}개 중 RNA-seq 코호트에도 측정된 유전자 {len(panel)}개')
print()


def test(d, genes, label):
    rows = []
    for g in genes:
        v = d.loc[g].dropna().values.astype(float)
        if len(v) < 8 or np.allclose(v, 0):
            continue
        rows.append((g, float(np.median(v)), 100 * float((v > 0).mean()), len(v),
                     stats.wilcoxon(v, alternative='two-sided').pvalue))
    r = pd.DataFrame(rows, columns=['gene', 'median_delta', 'pct_up', 'n', 'p']).set_index('gene')
    r['q'] = multipletests(r.p, method='fdr_bh')[1]
    print(f'{label}: {len(r)} genes · q<0.05 {int((r.q < 0.05).sum())} · q<0.10 {int((r.q < 0.10).sum())}')
    return r


RA = test(dA, panel, 'A GSE339484 RNA-seq 14쌍 · 패널로 한정')
RB = test(dB, panel, 'B GSE254874 NanoString 17쌍 · 같은 유전자')
print()

j = pd.DataFrame({'A_delta': RA.median_delta, 'A_q': RA.q,
                  'B_delta': RB.median_delta, 'B_q': RB.q}).dropna()
rho, pv = stats.spearmanr(j.A_delta, j.B_delta)
agree = ((j.A_delta > 0) == (j.B_delta > 0)).mean() * 100
print('=== 두 코호트 일치 ===')
print(f'  공통 {len(j)}개 · 변화량 Spearman rho={rho:.3f} (p={pv:.3g}) · 방향 일치 {agree:.1f}%')
both = j[(j.A_q < 0.05) & (j.B_q < 0.05)]
bothdir = both[(both.A_delta > 0) == (both.B_delta > 0)]
print(f'  양쪽 q<0.05 {len(both)}개 · 그중 방향까지 같은 것 {len(bothdir)}개')
print(f'  {sorted(bothdir.index)}')
j.sort_values('B_q').to_csv(f'{OUT}/replication.tsv', sep='\t')
print()

# B 의 상위 신호가 A 에서 같은 방향인가
topB = RB.sort_values('p').head(30).index
sub = j.loc[[g for g in topB if g in j.index]]
print(f'=== B 상위 30 중 A 에서 방향 일치 {int(((sub.A_delta>0)==(sub.B_delta>0)).sum())}/{len(sub)} ===')
print(sub.assign(같은방향=lambda x: np.where((x.A_delta > 0) == (x.B_delta > 0), 'O', 'X'))
      .head(18).to_string(float_format=lambda v: f'{v:.3f}'))
print()

# ---------------------------------------------------------------- 시간인가 생물학인가
soft = f'{D}/GSE254874_relapse.tsv'
if not os.path.exists(soft):
    import gzip, re
    t = gzip.decompress(urllib.request.urlopen(
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254874/soft/'
        'GSE254874_family.soft.gz', timeout=120).read()).decode('utf8', 'ignore')
    pid = re.findall(r'^!Sample_characteristics_ch1 = patientid:\s*(.*)$', t, re.M)
    rel = re.findall(r'^!Sample_characteristics_ch1 = relapse time \(d\):\s*(.*)$', t, re.M)
    pd.DataFrame({'patient': pid, 'relapse_d': rel}).drop_duplicates().to_csv(soft, sep='\t', index=False)
rt = pd.read_csv(soft, sep='\t')
rt['relapse_d'] = pd.to_numeric(rt.relapse_d, errors='coerce')
rt = rt.dropna().drop_duplicates('patient').set_index(rt.dropna().drop_duplicates('patient')
                                                      .patient.astype(str))
print('=== 재발까지 걸린 시간 ===')
print(f'  중앙 {rt.relapse_d.median():.0f}일 · 범위 {rt.relapse_d.min():.0f}–{rt.relapse_d.max():.0f}일 (n={len(rt)})')
common = [p for p in dB.columns if str(p) in rt.index]
if common:
    score = dB.loc[[g for g in RB[RB.q < 0.05].index], common].mean()
    days = rt.loc[[str(p) for p in common], 'relapse_d'].values
    rho2, pv2 = stats.spearmanr(score.values, days)
    print(f'  변화 크기 대 재발기간 Spearman rho={rho2:.3f} (p={pv2:.3g}, n={len(common)})')
    print('  → 유의하지 않으면 「시간이 흘러서」가 아니라 계통적 변화라는 뜻')
