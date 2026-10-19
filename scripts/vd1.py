#!/usr/bin/env python3
"""VitalDB: novelty and feasibility together, before any claim.

The candidate: intraoperative hypotension is defined everywhere by one universal
threshold (MAP < 65 mmHg). VitalDB holds four surgical specialties collected in
one hospital on one protocol, so the threshold can be tested for the same insult
across specialties without a batch difference -- the comparison single-specialty
studies cannot make.

Check the literature first, then check whether the registry actually carries what
the question needs.
"""
import urllib.request, urllib.parse, json, time, gzip, io
import pandas as pd, numpy as np

B = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def probe(label, q, k=4):
    u = B + 'esearch.fcgi?' + urllib.parse.urlencode(
        {'db': 'pubmed', 'term': q, 'retmax': k, 'retmode': 'json', 'sort': 'relevance'})
    r = json.load(urllib.request.urlopen(u, timeout=45))['esearchresult']
    print(f'### {label}: {int(r["count"])}건')
    if r['idlist']:
        u2 = B + 'esummary.fcgi?' + urllib.parse.urlencode(
            {'db': 'pubmed', 'id': ','.join(r['idlist']), 'retmode': 'json'})
        s = json.load(urllib.request.urlopen(u2, timeout=45))['result']
        for i in r['idlist']:
            print(f'    {s[i].get("pubdate","")[:4]} {s[i].get("source","")[:20]:22s} '
                  f'{s[i].get("title","")[:88]}')
    print()
    time.sleep(0.35)


probe('① 수술 중 저혈압 역치가 수술 종류에 따라 다른가',
      '("intraoperative hypotension") AND (threshold OR "MAP 65") AND '
      '("surgery type" OR "surgical specialty" OR "across procedures" OR heterogen*)')
probe('② 저혈압 역치 자체의 재검증',
      '("intraoperative hypotension") AND (threshold OR cut-off) AND (outcome) AND '
      '(nonlinear OR "change point" OR spline OR arbitrary)')
probe('③ VitalDB 를 쓴 연구',
      'VitalDB OR "vital recorder"')

# ---------------------------------------------------------------- feasibility
print('=' * 68)
print('=== VitalDB 가 실제로 무엇을 가지고 있나 ===')
raw = urllib.request.urlopen('https://api.vitaldb.net/cases', timeout=180).read()
if raw[:2] == b'\x1f\x8b':
    raw = gzip.decompress(raw)
d = pd.read_csv(io.BytesIO(raw), low_memory=False)
print(f'  증례 {len(d):,} · 변수 {len(d.columns)}')
print(f'  전문과: {dict(d.department.value_counts())}')

need = ['caseid', 'department', 'optype', 'age', 'sex', 'bmi', 'asa', 'emop',
        'death_inhosp', 'icu_days', 'los_postop', 'preop_cr', 'preop_hb', 'opend', 'opstart']
have = [c for c in need if c in d.columns]
print(f'  핵심 변수 확보: {have}')
print(f'  없는 것: {[c for c in need if c not in d.columns]}')

if 'death_inhosp' in d:
    print(f'\n  원내사망 {int(d.death_inhosp.sum())}건 ({100*d.death_inhosp.mean():.2f}%)')
    print('  전문과별 사망:')
    print(d.groupby('department').death_inhosp.agg(['size', 'sum', 'mean'])
          .rename(columns={'size': 'n', 'sum': '사망', 'mean': '비율'}).to_string())
if 'opstart' in d and 'opend' in d:
    dur = (d.opend - d.opstart) / 60
    print(f'\n  수술시간 중앙 {dur.median():.0f}분 (범위 {dur.min():.0f}-{dur.max():.0f})')
print(f'\n  트랙(파형) 목록 확인:')
tr = urllib.request.urlopen('https://api.vitaldb.net/trks', timeout=180).read()
if tr[:2] == b'\x1f\x8b':
    tr = gzip.decompress(tr)
t = pd.read_csv(io.BytesIO(tr), low_memory=False)
print(f'    트랙 레코드 {len(t):,} · 고유 신호 {t.tname.nunique()}')
key = [x for x in t.tname.unique() if any(k in x.upper() for k in
       ['ART_MBP', 'NIBP_MBP', 'ART_SBP', 'HR', 'SPO2'])]
print(f'    혈압·심박 관련: {key[:12]}')
cnt = t[t.tname.isin(key)].groupby('tname').caseid.nunique().sort_values(ascending=False)
print(f'    신호별 증례수:\n{cnt.head(8).to_string()}')
