#!/usr/bin/env python3
"""Source data for the expanded figure set (DK figure plan, 2026-10-14).

Starts from the deposited files, reproduces the 2026-09-15 pipeline (a1/a2/a3/gb2) and
asserts the confirmed numbers before computing anything new. New analyses requested in
FIG_PLAN_20261014.md, all kept inside the design rule (expression values are never
compared across cohorts; only within-patient differences, or statistics derived from
within-cohort data, cross between cohorts):

  PCA per cohort, sample QC, CAT curves, sign-flip and gene-label permutation nulls,
  patient bootstrap, elastic-net classifier with leave-one-pair-out, preranked GSEA
  (MSigDB Hallmark v2024.1), IDH-mutant exclusion, normalisation alternatives,
  pair-influence, panel genes removed by the cohort-A expression filter, module
  definition sensitivity.

Hyper-parameters are fixed in advance (no tuning). Seeds are fixed.
Writes data/derived/*.tsv and data/derived/summary.json.
"""
import json, os, re, gzip
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

H = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(H, '..')
OUT = os.path.join(P, 'data', 'derived'); os.makedirs(OUT, exist_ok=True)
S = {}                                   # summary numbers
rng = np.random.default_rng(20261014)

# ------------------------------------------------------------------ load (as a1)
rawA = pd.read_csv(f'{H}/GSE339484_counts_All_Samples.txt.gz', sep='\t', index_col=0)
logA_all = np.log2(rawA + 1)
keep = (logA_all > 1).sum(axis=1) >= logA_all.shape[1] * 0.5
logA = logA_all[keep]
assert logA.shape[0] == 23084
sym = json.load(open(f'{H}/out/ens2sym.json'))
def to_sym(df):
    d = df.copy(); d.index = [sym.get(g.split('.')[0], '') for g in d.index]
    return d[d.index != ''].groupby(level=0).mean()
A_expr = to_sym(logA)                                   # symbol x 28 samples
A_expr_all = to_sym(logA_all)                           # unfiltered, for S8
patA = lambda c: c.split('_')[0]; armA = lambda c: c.split('_')[1].lower()

b = pd.read_csv(f'{H}/GSE254874_processed.txt.gz', sep='\t', index_col=0)
def patB(c):
    c = c.strip(); return c.split()[0] if ' ' in c else c[:-1]
def armB(c):
    c = c.strip().lower()
    return c.split()[1] if ' ' in c else {'p': 'primary', 'r': 'recurrent'}.get(c[-1], c[-1])
B_expr = b

def pairs(df, pat, arm):
    pats = sorted({pat(c) for c in df.columns}); out = []
    for p_ in pats:
        pr = [c for c in df.columns if pat(c) == p_ and arm(c) == 'primary']
        rc = [c for c in df.columns if pat(c) == p_ and arm(c) == 'recurrent']
        if len(pr) == 1 and len(rc) == 1: out.append((p_, pr[0], rc[0]))
    return out
prA, prB = pairs(A_expr, patA, armA), pairs(B_expr, patB, armB)
assert len(prA) == 14 and len(prB) == 17
def deltas(df, pr): return pd.DataFrame({p_: df[r] - df[q] for p_, q, r in pr})
dA, dB = deltas(A_expr, prA), deltas(B_expr, prB)
shared = sorted(set(dA.index) & set(dB.index)); assert len(shared) == 699

def gene_tests(d, genes):
    rows = []
    for g in genes:
        v = d.loc[g].dropna().values.astype(float)
        if len(v) < 8 or np.allclose(v, 0): continue
        rows.append((g, float(np.median(v)), 100 * float((v > 0).mean()), len(v), stats.wilcoxon(v).pvalue))
    r = pd.DataFrame(rows, columns=['gene', 'median_delta', 'pct_up', 'n', 'p']).set_index('gene')
    r['q'] = multipletests(r.p, method='fdr_bh')[1]; return r
RA, RB = gene_tests(dA, shared), gene_tests(dB, shared)
j = pd.DataFrame({'A_delta': RA.median_delta, 'A_p': RA.p, 'A_q': RA.q, 'A_pct_up': RA.pct_up,
                  'B_delta': RB.median_delta, 'B_p': RB.p, 'B_q': RB.q, 'B_pct_up': RB.pct_up}).dropna()
ref = pd.read_csv(f'{H}/out/replication.tsv', sep='\t', index_col=0)
assert np.allclose(j.loc[ref.index, ['A_delta', 'B_delta', 'A_q', 'B_q']].values,
                   ref[['A_delta', 'B_delta', 'A_q', 'B_q']].values)
sig = j[j.B_q < 0.05].index; top30 = j.sort_values('B_q').head(30).index
def agree(a, bb, idx): return int(((a.loc[idx] > 0) == (bb.loc[idx] > 0)).sum())
assert agree(j.A_delta, j.B_delta, sig) == 64 and len(sig) == 67 and agree(j.A_delta, j.B_delta, top30) == 29
rho_obs = stats.spearmanr(j.A_delta, j.B_delta)[0]; assert round(rho_obs, 3) == 0.563
S['confirmed'] = dict(agree_sig='64/67', binom_p=stats.binomtest(64, 67).pvalue, agree_top30='29/30',
                      binom_p30=stats.binomtest(29, 30).pvalue, rho=rho_obs,
                      rho_p=stats.spearmanr(j.A_delta, j.B_delta)[1], agree_all=agree(j.A_delta, j.B_delta, j.index))
j.to_csv(f'{OUT}/TableS1_pergene_699.tsv', sep='\t')

# genome-wide cohort A (saved a1 result) with symbols and mean expression, for volcano/MA
gA = pd.read_csv(f'{H}/out/A_GSE339484_pergene.tsv', sep='\t', index_col=0)
gA['symbol'] = [sym.get(g.split('.')[0], '') for g in gA.index]
gA['mean_log2'] = logA.loc[gA.index].mean(axis=1)
gA.to_csv(f'{OUT}/A_genomewide.tsv', sep='\t')
S['A_genomewide'] = dict(n=len(gA), p_lt_001=int((gA.p < 0.01).sum()), expected=round(len(gA) * 0.01),
                         min_p=float(gA.p.min()), min_q=float(gA.q.min()))
gB = gene_tests(dB, list(dB.index)); gB['mean_log2'] = B_expr.loc[gB.index].mean(axis=1)
gB.to_csv(f'{OUT}/B_panel784.tsv', sep='\t')
S['B_panel784'] = dict(n=len(gB), q_lt_005=int((gB.q < 0.05).sum()))

# ------------------------------------------------------------------ metadata
soft = gzip.open('/tmp/claude-0/-workspace-icbmehs-c5-gbm/52e16d07-778b-458e-b957-5e9252622c4e/scratchpad/dl_geo_soft/GSE339484_family.soft.gz', 'rt').read()
metaA = {}
for blk in re.split(r'^\^SAMPLE = ', soft, flags=re.M)[1:]:
    t = re.findall(r'^!Sample_title = (.*)$', blk, re.M)[0]
    pid, grade, idh, _ = [x.strip() for x in t.split(',')]
    metaA[pid] = (grade, idh)
rel = pd.read_csv(f'{H}/GSE254874_relapse.tsv', sep='\t'); rel.index = rel.patient.astype(str)
meta = [dict(cohort='A', patient=p_, grade=metaA[p_][0], idh=metaA[p_][1], relapse_days='',
             primary_sample=q, recurrent_sample=r) for p_, q, r in prA]
meta += [dict(cohort='B', patient=p_, grade='Grade 4', idh='IDH wt', relapse_days=int(rel.loc[str(p_), 'relapse_d']),
              primary_sample=q, recurrent_sample=r) for p_, q, r in prB]
pd.DataFrame(meta).to_csv(f'{OUT}/TableS3_pairs.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ PCA + QC
def pca(X):                      # X genes x samples; centred per gene
    Z = (X.T - X.T.mean()).values
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    return pd.DataFrame(U[:, :3] * s[:3], index=X.columns, columns=['PC1', 'PC2', 'PC3']), s ** 2 / (s ** 2).sum()
for lab, X, pr, pat, arm in [('A', logA, prA, patA, armA), ('B', B_expr, prB, patB, armB)]:
    sc, ve = pca(X)
    sc['patient'] = [pat(c) for c in sc.index]; sc['arm'] = [arm(c) for c in sc.index]
    sc.to_csv(f'{OUT}/pca_{lab}.tsv', sep='\t')
    # how much of the variance is between patients: distance to own partner vs to others
    S[f'pca_{lab}'] = dict(var_explained=[float(x) for x in ve[:3]])
qcA = pd.DataFrame({'library_M': rawA.sum() / 1e6, 'detected_genes': (rawA > 0).sum(),
                    'expressed_genes_log2gt1': (logA_all > 1).sum()})
qcA['patient'] = [patA(c) for c in qcA.index]; qcA['arm'] = [armA(c) for c in qcA.index]
qcA.to_csv(f'{OUT}/qc_A.tsv', sep='\t')
qcB = pd.DataFrame({'median_log2cpm': B_expr.median(), 'iqr_log2cpm': B_expr.quantile(.75) - B_expr.quantile(.25)})
qcB['patient'] = [patB(c) for c in qcB.index]; qcB['arm'] = [armB(c) for c in qcB.index]
qcB.to_csv(f'{OUT}/qc_B.tsv', sep='\t')
S['qc'] = dict(A_lib_wilcoxon_p=float(stats.wilcoxon(qcA[qcA.arm == 'recurrent'].library_M.values -
                                                     qcA[qcA.arm == 'primary'].library_M.values).pvalue),
               A_lib_range=[float(qcA.library_M.min()), float(qcA.library_M.max())])
# sample-sample correlation within cohort (S1)
for lab, X in [('A', logA), ('B', B_expr)]:
    X.corr(method='spearman').to_csv(f'{OUT}/samplecorr_{lab}.tsv', sep='\t')

# ------------------------------------------------------------------ CAT curves
def signed(r): return np.sign(r.median_delta) * -np.log10(r.p)
sa, sb = signed(RA).loc[j.index], signed(RB).loc[j.index]
cat = []
for k in range(1, 301):
    upA, upB = set(sa.nlargest(k).index), set(sb.nlargest(k).index)
    dnA, dnB = set(sa.nsmallest(k).index), set(sb.nsmallest(k).index)
    cat.append((k, len(upA & upB) / k, len(dnA & dnB) / k, k / len(j)))
pd.DataFrame(cat, columns=['k', 'up', 'down', 'expected']).to_csv(f'{OUT}/cat.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ permutation nulls
DA = dA.loc[j.index].values; DB = dB.loc[j.index].values
bsig = j.B_delta.loc[sig].values > 0; isig = np.array([g in set(sig) for g in j.index])
NPERM = 10000
flip_agree, flip_rho, lab_agree, lab_rho = [], [], [], []
Bmed = j.B_delta.values
for _ in range(NPERM):
    s_ = rng.choice([-1, 1], size=DA.shape[1])
    am = np.median(DA * s_, axis=1)
    flip_agree.append(int(((am[isig] > 0) == bsig).sum())); flip_rho.append(stats.spearmanr(am, Bmed)[0])
    pm = rng.permutation(j.A_delta.values)
    lab_agree.append(int(((pm[isig] > 0) == bsig).sum())); lab_rho.append(stats.spearmanr(pm, Bmed)[0])
nulls = pd.DataFrame(dict(signflip_agree=flip_agree, signflip_rho=flip_rho, genelabel_agree=lab_agree, genelabel_rho=lab_rho))
nulls.to_csv(f'{OUT}/permutation_nulls.tsv', sep='\t', index=False)
pv = lambda x, o: (1 + int((np.asarray(x) >= o).sum())) / (len(x) + 1)
S['permutation'] = dict(n=NPERM, signflip_agree_p=pv(flip_agree, 64), signflip_rho_p=pv(flip_rho, rho_obs),
                        genelabel_agree_p=pv(lab_agree, 64), genelabel_rho_p=pv(lab_rho, rho_obs),
                        signflip_agree_max=int(max(flip_agree)), signflip_rho_max=float(max(flip_rho)),
                        signflip_agree_median=float(np.median(flip_agree)))

# ------------------------------------------------------------------ bootstrap (patients, within cohort)
NB = 2000; boot = []
itop = np.array([g in set(top30) for g in j.index]); btop = j.B_delta.values[itop] > 0
for _ in range(NB):
    ia = rng.integers(0, DA.shape[1], DA.shape[1]); ib = rng.integers(0, DB.shape[1], DB.shape[1])
    am, bm = np.median(DA[:, ia], axis=1), np.median(DB[:, ib], axis=1)
    boot.append((((am[isig] > 0) == (bm[isig] > 0)).mean(), ((am[itop] > 0) == (bm[itop] > 0)).mean(),
                 stats.spearmanr(am, bm)[0]))
boot = pd.DataFrame(boot, columns=['agree_sig', 'agree_top30', 'rho']); boot.to_csv(f'{OUT}/bootstrap.tsv', sep='\t', index=False)
S['bootstrap'] = {c: [float(np.percentile(boot[c], 2.5)), float(np.percentile(boot[c], 97.5))] for c in boot}

# ------------------------------------------------------------------ cross-cohort per-patient prediction (gb2) + full null
A699, B699 = dA.loc[j.index], dB.loc[j.index]
def cross(tr, te): d = tr.mean(axis=1).values; return np.array([stats.spearmanr(d, te[c].values)[0] for c in te.columns])
rAB, rBA = cross(A699, B699), cross(B699, A699)
assert round(np.median(rAB), 3) == 0.239 and round(np.median(rBA), 3) == 0.281
def cross_null(tr, te, n=2000):
    d = tr.mean(axis=1).values
    return np.array([np.median([stats.spearmanr(rng.permutation(d), te[c].values)[0] for c in te.columns]) for _ in range(n)])
nAB, nBA = cross_null(A699, B699), cross_null(B699, A699)
pd.DataFrame({'patient': list(B699.columns) + list(A699.columns), 'direction': ['A->B'] * 17 + ['B->A'] * 14,
              'rho': list(rAB) + list(rBA)}).to_csv(f'{OUT}/cross_patient.tsv', sep='\t', index=False)
pd.DataFrame({'null_median_AB': nAB, 'null_median_BA': nBA}).to_csv(f'{OUT}/cross_null_allpatients.tsv', sep='\t', index=False)
S['cross'] = dict(AB_median=float(np.median(rAB)), BA_median=float(np.median(rBA)),
                  AB_pos=int((rAB > 0).sum()), BA_pos=int((rBA > 0).sum()),
                  AB_null95=[float(np.percentile(nAB, 2.5)), float(np.percentile(nAB, 97.5))],
                  BA_null95=[float(np.percentile(nBA, 2.5)), float(np.percentile(nBA, 97.5))],
                  AB_null_p=pv(nAB, np.median(rAB)), BA_null_p=pv(nBA, np.median(rBA)))

# ------------------------------------------------------------------ elastic-net classifier, leave-one-pair-out
def model(): return LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=1.0, max_iter=50000)
def zfit(X): mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1; return mu, sd
def design(expr, pr, centre_patient=False):
    cols = [c for _, q, r in pr for c in (q, r)]
    X = expr.loc[j.index, cols].T.values.astype(float); y = np.array([0, 1] * len(pr)); g = np.repeat(np.arange(len(pr)), 2)
    if centre_patient:              # label-free: each patient's two samples centred on their own mean
        for k in range(len(pr)): X[g == k] -= X[g == k].mean(0)
    return X, y, g
def paired_acc(p, g): return float(np.mean([p[(g == k)][1] > p[(g == k)][0] for k in np.unique(g)]))
ml = {}; roc = []; lopo_rows = []
for variant, cp in [('sample', False), ('patient_centred', True)]:
    XA, yA, gA_ = design(A_expr, prA, cp); XB, yB, gB_ = design(B_expr, prB, cp)
    for lab, X, y, g, pr in [('A', XA, yA, gA_, prA), ('B', XB, yB, gB_, prB)]:
        pred = np.zeros(len(y))
        for k in np.unique(g):
            tr, te = g != k, g == k
            mu, sd = zfit(X[tr]); m = model().fit((X[tr] - mu) / sd, y[tr])
            pred[te] = m.predict_proba((X[te] - mu) / sd)[:, 1]
        auc = roc_auc_score(y, pred); ml[f'{variant}_LOPO_{lab}'] = dict(auc=auc, paired_acc=paired_acc(pred, g))
        fpr, tpr, _ = roc_curve(y, pred); roc += [(variant, f'LOPO {lab}', a, t) for a, t in zip(fpr, tpr)]
        lopo_rows += [(variant, lab, pr[k][0], pred[2 * k], pred[2 * k + 1]) for k in range(len(pr))]
    for lab, Xtr, ytr, Xte, yte, gte in [('A->B', XA, yA, XB, yB, gB_), ('B->A', XB, yB, XA, yA, gA_)]:
        mu, sd = zfit(Xtr); m = model().fit((Xtr - mu) / sd, ytr)
        mu2, sd2 = zfit(Xte)                         # test cohort standardised on itself, label-free
        pred = m.predict_proba((Xte - mu2) / sd2)[:, 1]
        auc = roc_auc_score(yte, pred)
        null = []
        for _ in range(2000):                        # swap labels within pairs
            sw = rng.integers(0, 2, len(np.unique(gte))).repeat(2).astype(bool)
            yy = np.where(sw, 1 - yte, yte); null.append(roc_auc_score(yy, pred))
        ml[f'{variant}_{lab}'] = dict(auc=auc, paired_acc=paired_acc(pred, gte), perm_p=pv(null, auc),
                                      n_nonzero=int((m.coef_ != 0).sum()))
        fpr, tpr, _ = roc_curve(yte, pred); roc += [(variant, lab, a, t) for a, t in zip(fpr, tpr)]
pd.DataFrame(roc, columns=['variant', 'curve', 'fpr', 'tpr']).to_csv(f'{OUT}/ml_roc.tsv', sep='\t', index=False)
pd.DataFrame(lopo_rows, columns=['variant', 'cohort', 'patient', 'p_primary', 'p_recurrent']).to_csv(f'{OUT}/ml_lopo_predictions.tsv', sep='\t', index=False)
S['ml'] = ml

# ------------------------------------------------------------------ modules
MODS = {'Myeloid/macrophage': ['CSF1R', 'SIRPA', 'MS4A6A', 'IRF8', 'DAB2', 'CMKLR1', 'AIF1', 'ITGAM', 'CD163', 'MRC1', 'TYROBP', 'FCGR3A', 'C1QA', 'C1QB'],
        'Stemness': ['SOX2', 'NOTCH1', 'OLIG2', 'NES', 'POU3F2', 'DLL3', 'ASCL1'],
        'Mesenchymal/ECM': ['VCAN', 'LOXL2', 'CDH11', 'FN1', 'COL1A1', 'SERPINE1', 'TGFB1', 'CD44']}
EXT = {'Myeloid/macrophage': MODS['Myeloid/macrophage'] + ['CD68', 'TREM2'],
       'Stemness': MODS['Stemness'] + ['PROM1', 'MYC', 'BMI1', 'FUT4', 'L1CAM'],
       'Mesenchymal/ECM': MODS['Mesenchymal/ECM'] + ['COL1A2', 'COL3A1', 'TIMP1', 'MMP2', 'MMP9', 'SPARC']}
def module_stats(d, gs):
    g = [x for x in gs if x in d.index]
    if len(g) < 3: return None
    v = d.loc[g].mean().values.astype(float)
    bm = [np.median(rng.choice(v, len(v))) for _ in range(2000)]
    return dict(n_genes=len(g), median=float(np.median(v)), ci=[float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))],
                pct_up=float(100 * (v > 0).mean()), p=float(stats.wilcoxon(v).pvalue), values=[float(x) for x in v], genes=g)
mod = {f'{m}|{lab}': module_stats(d, gs) for m, gs in MODS.items() for lab, d in [('A', dA), ('B', dB)]}
assert round(mod['Myeloid/macrophage|B']['median'], 3) == 1.042 and round(mod['Stemness|A']['median'], 3) == -0.994
assert round(mod['Mesenchymal/ECM|B']['median'], 3) == -0.711 and mod['Stemness|B'] is None
S['modules'] = {k: ({kk: vv for kk, vv in v.items() if kk != 'values'} if v else None) for k, v in mod.items()}
json.dump(mod, open(f'{OUT}/modules.json', 'w'))
# per-sample module scores (paired boxplots)
ms = []
for m, gs in MODS.items():
    for lab, X, pr in [('A', A_expr, prA), ('B', B_expr, prB)]:
        g = [x for x in gs if x in X.index]
        if len(g) < 3: continue
        for p_, q, r in pr: ms.append((m, lab, p_, float(X.loc[g, q].mean()), float(X.loc[g, r].mean())))
pd.DataFrame(ms, columns=['module', 'cohort', 'patient', 'primary', 'recurrent']).to_csv(f'{OUT}/module_sample_scores.tsv', sep='\t', index=False)
modrows = [dict(module=m, gene=g, in_A=g in dA.index, in_B=g in dB.index, core=g in MODS[m]) for m in EXT for g in EXT[m]]
pd.DataFrame(modrows).to_csv(f'{OUT}/TableS2_modules.tsv', sep='\t', index=False)
# S9: leave-one-gene-out and extended sets
s9 = []
for m in MODS:
    for lab, d in [('A', dA), ('B', dB)]:
        base = mod[f'{m}|{lab}']
        if base is None: continue
        for g in base['genes']:
            gs = [x for x in base['genes'] if x != g]; v = d.loc[gs].mean().values
            s9.append((m, lab, f'without {g}', len(gs), float(np.median(v)), float(stats.wilcoxon(v).pvalue)))
        e = module_stats(d, EXT[m])
        s9.append((m, lab, 'extended set', e['n_genes'], e['median'], e['p']))
pd.DataFrame(s9, columns=['module', 'cohort', 'variant', 'n_genes', 'median', 'p']).to_csv(f'{OUT}/S9_module_sensitivity.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ preranked GSEA (Hallmark)
gmt = {l.split('\t')[0]: l.strip().split('\t')[2:] for l in open(os.path.join(P, 'data', 'external', 'h.all.v2024.1.Hs.symbols.gmt'))}
def es_pos(pos, w, N):
    """Running-sum enrichment score from hit positions (weight p = 1), exact at the extremes."""
    pos = np.sort(pos); nh = len(pos); wh = w[pos]; cw = np.cumsum(wh) / wh.sum(); k = np.arange(nh)
    after = cw - (pos + 1 - (k + 1)) / (N - nh)
    before = np.concatenate([[0.0], cw[:-1]]) - (pos - k) / (N - nh)
    d = np.concatenate([after, before]); return d[np.argmax(np.abs(d))]
def gsea(score, nperm=2000, minsize=15):
    score = score.sort_values(ascending=False); genes = score.index.values; w = np.abs(score.values)
    N = len(genes); where = {g: i for i, g in enumerate(genes)}; rows = []
    for name, gs in gmt.items():
        pos = np.array([where[g] for g in set(gs) if g in where])
        if len(pos) < minsize: continue
        e = es_pos(pos, w, N)
        null = np.array([es_pos(rng.choice(N, len(pos), replace=False), w, N) for _ in range(nperm)])
        same = null[np.sign(null) == np.sign(e)]
        nes = e / np.abs(same).mean() if len(same) else np.nan
        p = (1 + (np.abs(same) >= abs(e)).sum()) / (1 + len(same))
        rows.append((name.replace('HALLMARK_', ''), len(pos), e, nes, p))
    r = pd.DataFrame(rows, columns=['set', 'size', 'ES', 'NES', 'p']).set_index('set')
    r['q'] = multipletests(r.p, method='fdr_bh')[1]; return r
gsym = gA[gA.symbol != ''].copy(); gsym['score'] = np.sign(gsym.median_delta) * -np.log10(gsym.p)
gsym = gsym.sort_values('p').drop_duplicates('symbol').set_index('symbol')
GA = gsea(gsym.score, minsize=15)
GB = gsea(np.sign(gB.median_delta) * -np.log10(gB.p), minsize=10)
GA.join(GB, lsuffix='_A', rsuffix='_B', how='outer').to_csv(f'{OUT}/gsea_hallmark.tsv', sep='\t')

# ------------------------------------------------------------------ S4: IDH-mutant exclusion
s4 = {}
idhmut = [p_ for p_, (g_, i_) in metaA.items() if 'mutant' in i_]
gr3wt = [p_ for p_, (g_, i_) in metaA.items() if 'mutant' not in i_ and g_ == 'Grade 3']
for name, drop in [('all 14 pairs', []), ('IDH-mutant excluded (11 pairs)', idhmut),
                   ('grade 4 IDH-wt only (10 pairs)', idhmut + gr3wt)]:
    dd = dA[[c for c in dA.columns if c not in drop]]
    r = gene_tests(dd, shared); jj = r.median_delta.loc[j.index]
    d699 = dd.loc[j.index]
    s4[name] = dict(n_pairs=dd.shape[1], agree_sig=agree(jj, j.B_delta, sig), agree_top30=agree(jj, j.B_delta, top30),
                    binom_p=stats.binomtest(agree(jj, j.B_delta, sig), 67).pvalue,
                    rho=float(stats.spearmanr(jj, j.B_delta)[0]), rho_p=float(stats.spearmanr(jj, j.B_delta)[1]),
                    A_q_lt_010=int((r.q < 0.10).sum()),
                    myeloid=module_stats(dd, MODS['Myeloid/macrophage']), stemness=module_stats(dd, MODS['Stemness']),
                    mes=module_stats(dd, MODS['Mesenchymal/ECM']),
                    cross_BA_median=float(np.median(cross(B699, d699))), cross_BA_pos=int((cross(B699, d699) > 0).sum()),
                    cross_AB_median=float(np.median(cross(d699, B699))), cross_AB_pos=int((cross(d699, B699) > 0).sum()))
    for k in ('myeloid', 'stemness', 'mes'):
        s4[name][k] = {kk: vv for kk, vv in s4[name][k].items() if kk not in ('values', 'genes')}
S['S4_idh'] = s4; S['S4_dropped'] = dict(idh_mutant=idhmut, grade3_wt=gr3wt)

# ------------------------------------------------------------------ S5: normalisation alternatives
def uq(raw):
    q = raw[raw.sum(axis=1) > 0].quantile(0.75); return np.log2(raw / q * q.mean() + 1)
def qnorm(X):
    r = X.rank(axis=0); m = np.sort(X.values, axis=0).mean(axis=1)
    return r.apply(lambda c: pd.Series(np.interp(c.values, np.arange(1, len(m) + 1), m), index=c.index))
s5 = []
altA = {'deposited, log2(x+1) (main)': A_expr, 'upper-quartile scaled': to_sym(uq(rawA)[keep]),
        'quantile normalised': to_sym(qnorm(logA))}
for nm, X in altA.items():
    dd = deltas(X, prA); r = gene_tests(dd, shared).median_delta.loc[j.index]
    s5.append(('A', nm, agree(r, j.B_delta, sig), 67, float(stats.spearmanr(r, j.B_delta)[0])))
Bc = B_expr - B_expr.median()
ddB = deltas(Bc, prB); rB2 = gene_tests(ddB, shared); sig2 = rB2[rB2.q < 0.05].index
s5.append(('B', 'per-sample median centred', agree(j.A_delta, rB2.median_delta, sig2), len(sig2),
           float(stats.spearmanr(j.A_delta, rB2.median_delta.loc[j.index])[0])))
pd.DataFrame(s5, columns=['cohort_varied', 'normalisation', 'agree', 'of', 'rho']).to_csv(f'{OUT}/S5_normalisation.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ S7: pair influence
s7 = []
for p_ in dA.columns:
    r = gene_tests(dA.drop(columns=p_), shared).median_delta.loc[j.index]
    s7.append(('A', p_, agree(r, j.B_delta, sig), 67, float(stats.spearmanr(r, j.B_delta)[0])))
for p_ in dB.columns:
    r = gene_tests(dB.drop(columns=p_), shared); sg = r[r.q < 0.05].index
    s7.append(('B', p_, agree(j.A_delta, r.median_delta, sg), len(sg), float(stats.spearmanr(j.A_delta, r.median_delta.loc[j.index])[0])))
pd.DataFrame(s7, columns=['cohort', 'left_out', 'agree', 'of', 'rho']).to_csv(f'{OUT}/S7_pair_influence.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ S8: panel genes removed by the A expression filter
outside = [g for g in dB.index if g not in shared]
dA_all = deltas(A_expr_all, prA)
rec = [g for g in outside if g in dA_all.index]
s8 = gB.loc[rec].copy(); s8['A_delta_unfiltered'] = dA_all.loc[rec].median(axis=1)
s8['A_mean_log2'] = A_expr_all.loc[rec].mean(axis=1)
s8.to_csv(f'{OUT}/S8_filtered_panel_genes.tsv', sep='\t')
sg8 = s8[s8.q < 0.05]
S['S8'] = dict(panel_not_shared=len(outside), recoverable_unfiltered=len(rec), not_in_A_at_all=len(outside) - len(rec),
               B_sig_among_recovered=len(sg8), agree=int(((sg8.median_delta > 0) == (sg8.A_delta_unfiltered > 0)).sum()))

json.dump(S, open(f'{OUT}/summary.json', 'w'), indent=1, default=float)
print(json.dumps({k: S[k] for k in ['permutation', 'bootstrap', 'cross', 'ml', 'S8', 'qc']}, indent=1, default=float))
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ('n_pairs', 'agree_sig', 'binom_p', 'rho', 'agree_top30', 'cross_AB_pos', 'cross_BA_pos')} for k, v in S['S4_idh'].items()}, indent=1, default=float))
