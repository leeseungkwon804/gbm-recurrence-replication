#!/usr/bin/env python3
"""The decisive test.

The Nature Communications paper reports two pathology-specific findings:
neuroinflammation in TLE-HS, and oligodendrocyte dysfunction in the
mTORopathies. Both are group contrasts between an adult temporal-lobe disease
and paediatric frontal-lobe diseases. Do they survive adjustment for age and
region?

Gene sets are canonical cell-type and process markers, fixed here before the
test and not chosen from the data.
"""
import gzip, re, urllib.request
import numpy as np, pandas as pd
from scipy import stats

D = '/workspace/_scout_ksh/ep'

NEUROINFLAM = ['AIF1', 'CD68', 'ITGAM', 'TYROBP', 'C1QA', 'C1QB', 'C1QC', 'C3', 'CD14',
               'FCGR3A', 'TLR2', 'TLR4', 'IL1B', 'TNF', 'CCL2', 'CXCL8', 'GFAP',
               'S100B', 'VIM', 'CD44', 'TREM2', 'CSF1R', 'P2RY12', 'CX3CR1']
OLIGO = ['MBP', 'PLP1', 'MOG', 'MAG', 'MOBP', 'CNP', 'OLIG1', 'OLIG2', 'SOX10',
         'CLDN11', 'MAL', 'UGT8', 'ERMN', 'ASPA', 'OPALIN', 'GJB1', 'CA2', 'TF']

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
m['ctrl'] = m.dz.str.contains('Control')

e = pd.read_csv(f'{D}/GSE256068_processed_data.csv.gz', sep=r'\s+', index_col=0,
                engine='python')
# the file keeps literal quote characters around every field
e.index = e.index.astype(str).str.strip('"')
e.columns = m.gsm.values

# ensembl -> symbol
import json, os, urllib.parse
cache = f'{D}/sym.json'
if os.path.exists(cache):
    s2 = json.load(open(cache))
else:
    s2 = {}
    ids = [g.split('.')[0] for g in e.index]
    for i in range(0, len(ids), 900):
        body = urllib.parse.urlencode({'ids': ','.join(ids[i:i + 900]),
                                       'fields': 'symbol', 'species': 'human'}).encode()
        try:
            for h in json.load(urllib.request.urlopen('https://mygene.info/v3/gene', body, timeout=120)):
                if h.get('symbol'):
                    s2[h['query']] = h['symbol']
        except Exception as ex:
            print('mygene 실패', ex); break
    json.dump(s2, open(cache, 'w'))
e = e.copy()
e['sym'] = [s2.get(g.split('.')[0], '') for g in e.index]
E = e[e.sym != ''].groupby('sym').mean()
print(f'유전자 {E.shape[0]:,} · 샘플 {E.shape[1]}')

dis = m[~m.ctrl].reset_index(drop=True)
dis['grp'] = np.where(dis.dz == 'TLE-HS', 'TLE-HS', 'mTORopathy')
print(f'TLE-HS {int((dis.grp=="TLE-HS").sum())} · mTOR병증 {int((dis.grp=="mTORopathy").sum())}')
print()

def score(genes):
    g = [x for x in genes if x in E.index]
    z = E.loc[g, dis.gsm.values]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0).values, len(g)


for nm, gs in [('신경염증 (NC: TLE-HS 특이)', NEUROINFLAM),
               ('희소돌기아교 (NC: mTOR 특이)', OLIGO)]:
    y, ng = score(gs)
    a = y[(dis.grp == 'TLE-HS').values]
    b = y[(dis.grp == 'mTORopathy').values]
    U, p_raw = stats.mannwhitneyu(a, b)
    # adjusted: regress out age + region, then retest on residuals
    reg = pd.get_dummies(dis.region, drop_first=True).astype(float).values
    age = dis.age.fillna(dis.age.median()).values.reshape(-1, 1)
    X = np.hstack([np.ones((len(dis), 1)), age, reg])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    p_adj = stats.mannwhitneyu(r[(dis.grp == 'TLE-HS').values],
                               r[(dis.grp == 'mTORopathy').values]).pvalue
    d_raw = np.median(a) - np.median(b)
    d_adj = (np.median(r[(dis.grp == 'TLE-HS').values])
             - np.median(r[(dis.grp == 'mTORopathy').values]))
    print(f'### {nm}   (유전자 {ng}개)')
    print(f'    보정 전: 차이 {d_raw:+.3f}  p = {p_raw:.3g}')
    print(f'    나이·부위 보정 후: 차이 {d_adj:+.3f}  p = {p_adj:.3g}')
    print(f'    → {"살아남음" if p_adj < 0.05 else "⚠ 사라짐"}')
    print()

# 측두엽·나이맞춤 부분집합에서도
w = dis[(dis.region == 'Temporal') & (dis.age.between(10, 40))]
print(f'=== 측두엽 + 10-40세 부분집합 (n={len(w)}: '
      f'{dict(w.grp.value_counts())}) ===')
for nm, gs in [('신경염증', NEUROINFLAM), ('희소돌기아교', OLIGO)]:
    g = [x for x in gs if x in E.index]
    z = E.loc[g, w.gsm.values]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
    y = z.mean(axis=0).values
    a = y[(w.grp == 'TLE-HS').values]; b = y[(w.grp == 'mTORopathy').values]
    if len(b) >= 3:
        p = stats.mannwhitneyu(a, b).pvalue
        print(f'  {nm}: 차이 {np.median(a)-np.median(b):+.3f}  p = {p:.3g}  '
              f'(n {len(a)} 대 {len(b)})')
