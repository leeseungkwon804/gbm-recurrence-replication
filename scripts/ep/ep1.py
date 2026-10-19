#!/usr/bin/env python3
"""Are pathology, age and brain region separable in GSE256068?

The Nature Communications analysis compares TLE-HS against the mTORopathies and
reports pathology-specific mechanisms. But TLE-HS has a median age of 37 and is
entirely temporal, while TSC and the FCDs are paediatric and mostly frontal. This
asks how much of the between-group variance each factor can claim once the others
are in the model, and whether an age- and region-matched comparison survives.
"""
import gzip, re, collections
import numpy as np, pandas as pd
import urllib.request

D = '/workspace/_scout_ksh/ep'

# ---------------------------------------------------------------- metadata
t = gzip.decompress(open(f'{D}/soft.gz', 'rb').read()).decode('utf8', 'ignore') \
    if False else gzip.decompress(urllib.request.urlopen(
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256068/soft/'
        'GSE256068_family.soft.gz', timeout=180).read()).decode('utf8', 'ignore')
rec = []
for b in t.split('^SAMPLE = ')[1:]:
    gsm = b.split('\n')[0].strip()
    ch = {}
    for c in re.findall(r'^!Sample_characteristics_ch1 = (.*)$', b, re.M):
        if ':' in c:
            k, v = c.split(':', 1)
            ch[k.strip().lower()] = v.strip()
    ti = re.search(r'^!Sample_title = (.*)$', b, re.M)
    rec.append(dict(gsm=gsm, title=ti.group(1) if ti else '', **ch))
m = pd.DataFrame(rec)


def num(x):
    try:
        return float(re.sub(r'[^0-9.\-]', '', str(x)))
    except Exception:
        return np.nan


m['freq'] = m['seisure frequency'].map(num)
m['age'] = m['age'].map(num)
m = m.rename(columns={'disease state': 'dz', 'tissue': 'region'})
m['is_ctrl'] = m.dz.str.contains('Control')
print(f'메타 {len(m)}행 · 컬럼 {list(m.columns)}')

# ---------------------------------------------------------------- expression
e = pd.read_csv(f'{D}/GSE256068_processed_data.csv.gz', sep=r'\s+', index_col=0,
                engine='python')
print(f'발현 {e.shape[0]:,} 유전자 × {e.shape[1]} 샘플')
if e.shape[1] == len(m):
    e.columns = m.gsm.values          # order assumed to follow the SOFT record order
print(f'샘플 수 일치: {e.shape[1] == len(m)}')

keep = (e > 1).sum(axis=1) >= e.shape[1] * 0.3
E = e[keep]
print(f'발현 검출 유전자 {E.shape[0]:,}')
print()

# ---------------------------------------------------------------- variance decomposition
def rss(Y, X):
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    r = Y - X @ b
    return (r ** 2).sum(axis=0)


sub = m[~m.is_ctrl].reset_index(drop=True)
Y = E[sub.gsm.values].T.values.astype(float)          # samples x genes
Y = Y - Y.mean(axis=0)
n = len(sub)
one = np.ones((n, 1))


def dummies(s):
    lv = sorted(s.unique())[1:]
    return np.column_stack([(s == l).astype(float) for l in lv]) if lv else np.zeros((n, 0))


blocks = {
    'pathology': dummies(sub.dz),
    'age': sub.age.fillna(sub.age.median()).values.reshape(-1, 1),
    'region': dummies(sub.region),
    'sex': dummies(sub.sex),
}
full = np.hstack([one] + [blocks[k] for k in blocks])
tss = (Y ** 2).sum(axis=0)
r_full = rss(Y, full)
print('=== 분산 분해 (질환 샘플 %d개, 유전자당 부분 R², 중앙값) ===' % n)
print(f'{"요인":<12}{"단독 R²":>10}{"고유 R²(나머지 통제)":>22}')
for k in blocks:
    solo = np.hstack([one, blocks[k]])
    r_solo = rss(Y, solo)
    red = np.hstack([one] + [blocks[j] for j in blocks if j != k])
    r_red = rss(Y, red)
    print(f'{k:<12}{np.median(1 - r_solo / tss):>10.3f}'
          f'{np.median((r_red - r_full) / tss):>22.3f}')
print(f'{"전체 모형":<12}{np.median(1 - r_full / tss):>10.3f}')
print()

# ---------------------------------------------------------------- matched subset
print('=== 나이·부위를 맞추면 몇 개나 남나 ===')
for lo, hi in [(10, 40), (15, 35), (2, 46)]:
    w = sub[(sub.age >= lo) & (sub.age <= hi)]
    tab = pd.crosstab(w.dz, w.region)
    print(f'  나이 {lo}-{hi}세: n={len(w)}')
    print(tab.to_string().replace('\n', '\n     '))
    tw = w[w.region == 'Temporal']
    print(f'     └ 측두엽만: {dict(tw.dz.value_counts())}')
    print()
