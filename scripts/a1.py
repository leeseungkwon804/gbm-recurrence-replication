#!/usr/bin/env python3
"""First pass for the Kim Seunghyun project: does glioblastoma change between the
first resection and re-operation, within the same patient?

Design rule taken from C2: the two cohorts sit on different platforms (RNA-seq
and an 784-gene targeted panel), so nothing is pooled at the level of expression.
Each patient contributes a within-person difference computed inside its own
dataset, and only those differences are compared between cohorts.

No gene set is imposed yet. The first question is the unbiased one: how much
changes at all, and do two independent cohorts agree about what changes?
"""
import gzip, os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

D = '/workspace/_scout_ksh'
OUT = f'{D}/out'
os.makedirs(OUT, exist_ok=True)


def paired(df, pat_of, arm_of):
    """genes x patients matrix of (recurrent - primary), inside one cohort."""
    pats = sorted({pat_of(c) for c in df.columns})
    cols = {}
    for p in pats:
        pri = [c for c in df.columns if pat_of(c) == p and arm_of(c) == 'primary']
        rec = [c for c in df.columns if pat_of(c) == p and arm_of(c) == 'recurrent']
        if len(pri) == 1 and len(rec) == 1:
            cols[p] = df[rec[0]] - df[pri[0]]
    return pd.DataFrame(cols)


# ---------------------------------------------------------------- cohort A: RNA-seq
a = pd.read_csv(f'{D}/GSE339484_counts_All_Samples.txt.gz', sep='\t', index_col=0)
# values are normalised counts; put them on a log scale like the panel cohort
a = np.log2(a + 1)
# drop genes that are essentially absent
keep = (a > 1).sum(axis=1) >= a.shape[1] * 0.5
a = a[keep]
dA = paired(a, lambda c: c.split('_')[0], lambda c: c.split('_')[1].lower())
print(f'A  GSE339484  {a.shape[0]:,} genes expressed · {dA.shape[1]} paired patients')

# ---------------------------------------------------------------- cohort B: 784-gene panel
b = pd.read_csv(f'{D}/GSE254874_processed.txt.gz', sep='\t', index_col=0)
# the last three pairs in this series use '599p'/'599r' instead of '599 primary'
def b_pat(c):
    c = c.strip()
    return c.split()[0] if ' ' in c else c[:-1]


def b_arm(c):
    c = c.strip().lower()
    if ' ' in c:
        return c.split()[1]
    return {'p': 'primary', 'r': 'recurrent'}.get(c[-1], c[-1])


dB = paired(b, b_pat, b_arm)
print(f'B  GSE254874  {b.shape[0]:,} panel genes · {dB.shape[1]} paired patients')
print()


def test(d, label):
    rows = []
    for g in d.index:
        v = d.loc[g].dropna().values
        if len(v) < 8 or np.allclose(v, 0):
            continue
        try:
            p = stats.wilcoxon(v, alternative='two-sided').pvalue
        except ValueError:
            continue
        rows.append((g, np.median(v), 100 * (v > 0).mean(), len(v), p))
    r = pd.DataFrame(rows, columns=['gene', 'median_delta', 'pct_up', 'n', 'p']).set_index('gene')
    r['q'] = multipletests(r.p, method='fdr_bh')[1]
    r = r.sort_values('p')
    print(f'{label}: {len(r):,} genes tested · q<0.05 {int((r.q < 0.05).sum())} '
          f'· q<0.10 {int((r.q < 0.10).sum())}')
    up = r[(r.q < 0.05) & (r.median_delta > 0)]
    dn = r[(r.q < 0.05) & (r.median_delta < 0)]
    print(f'   올라간 것 {len(up)} · 내려간 것 {len(dn)}')
    return r


rA = test(dA, 'A GSE339484 (RNA-seq, 14쌍)')
rB = test(dB, 'B GSE254874 (panel, 17쌍)')
print()

# ---------------------------------------------------------------- do the two agree?
def symbols(ens):
    """Ensembl gene id -> symbol, in batches, from mygene.info."""
    import urllib.request, urllib.parse, json
    out = {}
    ens = [e.split('.')[0] for e in ens]
    for i in range(0, len(ens), 500):
        chunk = ens[i:i + 500]
        data = urllib.parse.urlencode({'ids': ','.join(chunk), 'fields': 'symbol',
                                       'species': 'human'}).encode()
        try:
            r = json.load(urllib.request.urlopen('https://mygene.info/v3/gene', data, timeout=90))
            for h in r:
                if h.get('symbol'):
                    out[h['query']] = h['symbol']
        except Exception as e:
            print('   symbol 조회 실패:', e)
            break
    return out


sigA = rA[rA.q < 0.05].index.tolist()
mapA = symbols(sigA[:2000])
rA['symbol'] = [mapA.get(g.split('.')[0], '') for g in rA.index]
shared = sorted(set(rA.loc[rA.q < 0.05, 'symbol']) & set(rB.index))
both = sorted(set(rA.symbol) & set(rB.index))
print(f'=== 두 코호트 겹침 ===')
print(f'   패널 784개 중 RNA-seq 에도 있는 유전자 {len(both)}개')
if both:
    j = pd.DataFrame({'A_delta': [rA.loc[rA.symbol == s, 'median_delta'].iloc[0] for s in both],
                      'A_q': [rA.loc[rA.symbol == s, 'q'].iloc[0] for s in both],
                      'B_delta': rB.loc[both, 'median_delta'].values,
                      'B_q': rB.loc[both, 'q'].values}, index=both)
    ok = j.dropna()
    rho, pv = stats.spearmanr(ok.A_delta, ok.B_delta)
    agree = ((ok.A_delta > 0) == (ok.B_delta > 0)).mean() * 100
    print(f'   변화량 상관 Spearman rho={rho:.3f} (p={pv:.3g}) · 방향 일치 {agree:.1f}%')
    bs = ok[(ok.A_q < 0.05) & (ok.B_q < 0.05)]
    print(f'   양쪽 모두 q<0.05 인 유전자 {len(bs)}개: {list(bs.index)[:20]}')
    j.to_csv(f'{OUT}/cohort_agreement.tsv', sep='\t')
print()

rA.to_csv(f'{OUT}/A_GSE339484_pergene.tsv', sep='\t')
rB.to_csv(f'{OUT}/B_GSE254874_pergene.tsv', sep='\t')
dA.to_csv(f'{OUT}/A_deltas.tsv', sep='\t')
dB.to_csv(f'{OUT}/B_deltas.tsv', sep='\t')

print('=== B 코호트 상위 15 (패널이라 해석이 쉽다) ===')
print(rB.head(15).to_string(float_format=lambda x: f'{x:.4g}'))
print()
print('=== 환자별 전체 변화량 (중앙 |Δ|) — 사람마다 얼마나 다른가 ===')
for lab, d in [('A', dA), ('B', dB)]:
    m = d.abs().median(axis=0).sort_values()
    print(f'  {lab}: 최소 {m.iloc[0]:.3f} · 중앙 {m.median():.3f} · 최대 {m.iloc[-1]:.3f}  (n={len(m)})')
