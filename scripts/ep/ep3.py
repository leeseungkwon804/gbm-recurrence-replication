#!/usr/bin/env python3
"""Does seizure burden leave a molecular trace, and is it the same trace in
different epilepsy pathologies?

Discovery inside TLE-HS only: 64 samples, every one temporal lobe, so pathology
and region are constants and cannot confound. Age and sex are adjusted. Seizure
frequency is heavily skewed so it enters as log.

Then the signature is tested in each of the other three pathologies separately,
which are different patients, different lobes and different ages -- a real
replication rather than a re-fit.

The Nature Communications paper that produced this dataset compared pathology
groups and never used the seizure-frequency field.
"""
import gzip, re, json, os, urllib.request
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

D = '/workspace/_scout_ksh/ep'

t = gzip.decompress(urllib.request.urlopen(
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256068/soft/'
    'GSE256068_family.soft.gz', timeout=180).read()).decode('utf8', 'ignore')
rec = []
for b in t.split('^SAMPLE = ')[1:]:
    ch = {}
    for c in re.findall(r'^!Sample_characteristics_ch1 = (.*)$', b, re.M):
        if ':' in c:
            k, v = c.split(':', 1)
            ch[k.strip().lower()] = v.strip()
    rec.append(dict(gsm=b.split('\n')[0].strip(), **ch))
m = pd.DataFrame(rec).rename(columns={'disease state': 'dz', 'tissue': 'region'})


def num(x):
    try:
        return float(re.sub(r'[^0-9.\-]', '', str(x)))
    except Exception:
        return np.nan


m['age'] = m.age.map(num)
m['freq'] = m['seisure frequency'].map(num)
m['male'] = (m.sex == 'Male').astype(float)

e = pd.read_csv(f'{D}/GSE256068_processed_data.csv.gz', sep=r'\s+', index_col=0,
                engine='python')
e.index = e.index.astype(str).str.strip('"')
e.columns = m.gsm.values
s2 = json.load(open(f'{D}/sym.json'))
e = e.copy()
e['sym'] = [s2.get(g.split('.')[0], '') for g in e.index]
E = e[e.sym != ''].groupby('sym').mean()
E = E[(E > 1).sum(axis=1) >= E.shape[1] * 0.3]
print(f'유전자 {E.shape[0]:,} · 샘플 {E.shape[1]}')

# ---------------------------------------------------------------- discovery: TLE-HS
d = m[(m.dz == 'TLE-HS') & m.freq.notna()].reset_index(drop=True)
print(f'\n=== 발굴 코호트 TLE-HS n={len(d)} · 부위 {set(d.region)} ===')
print(f'    발작빈도 중앙 {d.freq.median():.0f} (범위 {d.freq.min():.0f}-{d.freq.max():.0f})'
      f' · 나이 중앙 {d.age.median():.0f}')
rho_af, p_af = stats.spearmanr(d.freq, d.age)
print(f'    이 군 안에서 빈도-나이 상관 rho={rho_af:.3f} (p={p_af:.3g}) '
      f'→ {"교란 주의" if p_af < 0.05 else "독립적"}')

Y = E[d.gsm.values].T.values.astype(float)
lf = np.log10(d.freq.values + 1)
X = np.column_stack([np.ones(len(d)), d.age.fillna(d.age.median()).values, d.male.values])


def resid(v, X):
    b, *_ = np.linalg.lstsq(X, v, rcond=None)
    return v - X @ b


Yr = resid(Y, X)
lfr = resid(lf.reshape(-1, 1), X).ravel()
rows = []
for j, g in enumerate(E.index):
    r, p = stats.spearmanr(Yr[:, j], lfr)
    rows.append((g, r, p))
res = pd.DataFrame(rows, columns=['gene', 'rho', 'p']).set_index('gene')
res['q'] = multipletests(res.p, method='fdr_bh')[1]
res = res.sort_values('p')
sig = res[res.q < 0.05]
print(f'    → 발작빈도와 연관된 유전자: q<0.05 {len(sig)}개 · q<0.10 {int((res.q<0.10).sum())}'
      f' · p<0.001 {int((res.p<0.001).sum())} (기대 {0.001*len(res):.0f})')
if len(sig):
    print(f'    상위: {list(sig.head(12).index)}')
    print(f'    양의 상관 {int((sig.rho>0).sum())} · 음의 상관 {int((sig.rho<0).sum())}')
res.to_csv(f'{D}/tle_freq_assoc.tsv', sep='\t')

# ---------------------------------------------------------------- replication
if len(sig) >= 10:
    up = list(sig[sig.rho > 0].index)
    dn = list(sig[sig.rho < 0].index)
    print('\n=== 다른 병리에서 같은 신호가 나오나 (독립 환자·다른 엽) ===')
    for dz in ['FCD2b', 'TSC', 'FCD2a']:
        w = m[(m.dz == dz) & m.freq.notna()].reset_index(drop=True)
        if len(w) < 10:
            print(f'  {dz}: n={len(w)} 너무 작음'); continue
        Z = E[w.gsm.values]
        z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1) + 1e-9, axis=0)
        score = (z.loc[[g for g in up if g in z.index]].mean(axis=0)
                 - z.loc[[g for g in dn if g in z.index]].mean(axis=0)).values
        Xw = np.column_stack([np.ones(len(w)), w.age.fillna(w.age.median()).values,
                              w.male.values])
        sr = resid(score.reshape(-1, 1), Xw).ravel()
        fr = resid(np.log10(w.freq.values + 1).reshape(-1, 1), Xw).ravel()
        r, p = stats.spearmanr(sr, fr)
        print(f'  {dz:7s} n={len(w):<3} 부위 {dict(w.region.value_counts())} → '
              f'점수-빈도 rho={r:+.3f} p={p:.3g} {"✓" if p < 0.05 else ""}')
