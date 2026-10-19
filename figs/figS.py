#!/usr/bin/env python3
"""Supplementary Figures S1-S10 (FIG_PLAN_20261014; layout per HQ_c5_FIG_FIX_20261014).
Output: figs/supp/FigureS*.pdf/png. Same grid, letters and font as the main figures.
S9 deviates from the plan: the 85 panel genes outside the shared 699 have no cohort-A
counterpart at all (0 recoverable even without the expression filter), so concordance
outside the shared set cannot be computed. S9 instead shows the cohort-A signal outside
the panel (genome-wide p-value distribution and Hallmark agreement)."""
import json, os
import numpy as np, pandas as pd
from scipy import stats
from style import *

SUP = os.path.join(HERE, 'supp'); os.makedirs(SUP, exist_ok=True)
S = json.load(open(f'{DER}/summary.json'))
j = pd.read_csv(f'{DER}/TableS1_pergene_699.tsv', sep='\t', index_col=0)
H = os.path.join(HERE, '..', 'scripts')
sym = json.load(open(f'{H}/out/ens2sym.json'))
dA = pd.read_csv(f'{H}/out/A_deltas.tsv', sep='\t', index_col=0)
dA.index = [sym.get(g.split('.')[0], '') for g in dA.index]; dA = dA[dA.index != ''].groupby(level=0).mean()
dB = pd.read_csv(f'{H}/out/B_deltas.tsv', sep='\t', index_col=0)
meta = pd.read_csv(f'{DER}/TableS3_pairs.tsv', sep='\t')
sig = j[j.B_q < 0.05]
OKALL = True
def done(L, name):
    global OKALL
    OKALL &= L.check(name); L.save(name.replace('Figure ', 'Figure'), SUP)

# ---------------------------------------------------------------- S1 QC
rawA = pd.read_csv(f'{H}/GSE339484_counts_All_Samples.txt.gz', sep='\t', index_col=0)
logA = np.log2(rawA + 1); logA = logA[(logA > 1).sum(axis=1) >= 14]
b = pd.read_csv(f'{H}/GSE254874_processed.txt.gz', sep='\t', index_col=0)
L = Layout([dict(h=1.55, top=0.26, bottom=0.42), dict(h=1.55, top=0.26, bottom=0.42), dict(h=2.75, top=0.26, bottom=0.30, right=0.62)])
for r, (X, lab, letter) in enumerate([(logA, 'A', 'a'), (b, 'B', 'b')]):
    ax = L.ax(r, 0, 12, letter=letter, title=f'Cohort {lab}: per-sample distribution ({X.shape[0]:,} genes)')
    ax.boxplot([X[c].values for c in X.columns], widths=0.6, showfliers=False, medianprops=dict(color=COH[lab], lw=1),
               boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    ax.set_xticks(range(1, X.shape[1] + 1), [c.replace('Patient', 'P').replace('_Primary', 'p').replace('_Recurrent', 'r').replace(' primary', 'p').replace(' recurrent', 'r').strip() for c in X.columns], rotation=90, fontsize=4.8)
    ax.set_ylabel('log$_2$ expression' if lab == 'A' else 'log$_2$ CPM')
for k, (lab, letter) in enumerate([('A', 'c'), ('B', 'd')]):
    C = pd.read_csv(f'{DER}/samplecorr_{lab}.tsv', sep='\t', index_col=0)
    ax = L.ax(2, 6 * k, 6, letter=letter, title=f'Sample–sample correlation, cohort {lab}', left=0.30)
    im = ax.imshow(C.values, cmap='Greys', vmin=C.values.min(), vmax=1, interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('Samples in deposited order (pairs adjacent)', fontsize=6)
    x, y, w, h = L.box(2, 6 * k, 6, left=0.30)
    cax = L.fig.add_axes([(x + w + 0.08) / L.W, y / L.H, 0.07 / L.W, h / L.H]); cb = L.fig.colorbar(im, cax=cax)
    cb.set_label('Spearman r', fontsize=6); cb.ax.tick_params(labelsize=5.5)
done(L, 'Figure S1')

# ---------------------------------------------------------------- S2 covariates
pcaA = pd.read_csv(f'{DER}/pca_A.tsv', sep='\t', index_col=0); pcaB = pd.read_csv(f'{DER}/pca_B.tsv', sep='\t', index_col=0)
qcA = pd.read_csv(f'{DER}/qc_A.tsv', sep='\t', index_col=0)
mA = meta[meta.cohort == 'A'].set_index('patient'); mB = meta[meta.cohort == 'B'].set_index('patient')
grp = {p_: ('IDH1-mutant' if 'mutant' in mA.loc[p_, 'idh'] else ('grade 3 IDH1-wt' if mA.loc[p_, 'grade'] == 'Grade 3' else 'grade 4 IDH1-wt')) for p_ in mA.index}
cols = {'grade 4 IDH1-wt': COH['A'], 'grade 3 IDH1-wt': '#e08214', 'IDH1-mutant': INK}
L = Layout([dict(h=2.25, top=0.30, bottom=0.62), dict(h=2.15, top=0.30, bottom=0.42, right=0.45)])
ax = L.ax(0, 0, 6, letter='a', title='Cohort A PCA by grade and IDH status')
for g_, c in cols.items():
    s_ = pcaA[pcaA.patient.map(grp) == g_]
    ax.scatter(s_.PC1, s_.PC2, s=12, facecolor=np.where(s_.arm == 'recurrent', c, 'white'), edgecolor=c, lw=0.8)
ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
L.legend_row(0, [plt.Line2D([], [], ls='none', marker='o', ms=4, mfc=c, mec=c) for c in cols.values()], list(cols), ncol=3, c0=0, span=6, y_off=0.0, fontsize=5.8)
r_, p_ = stats.spearmanr(qcA.library_M, pcaA.loc[qcA.index, 'PC1'])
ax = L.ax(0, 6, 6, letter='b', title='Cohort A: PC1 vs column sum', note=f'ρ = {r_:.2f}, p = {sci(p_)}')
ax.scatter(qcA.library_M, pcaA.loc[qcA.index, 'PC1'], s=12, facecolor=np.where(qcA.arm == 'recurrent', COH['A'], 'white'), edgecolor=COH['A'], lw=0.8)
ax.set_xlabel('Column sum of deposited values (M)'); ax.set_ylabel('PC1')
ax = L.ax(1, 0, 6, letter='c', title='Cohort B PCA by relapse interval')
rel = pcaB.patient.astype(str).map(mB.relapse_days.rename(index=str))
sc = ax.scatter(pcaB.PC1, pcaB.PC2, c=rel, cmap='Greys', vmin=0, s=12, edgecolor=COH['B'], lw=0.6)
ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
x, y, w, h = L.box(1, 0, 6); cax = L.fig.add_axes([(x + w + 0.08) / L.W, y / L.H, 0.07 / L.W, h / L.H])
cb = L.fig.colorbar(sc, cax=cax); cb.set_label('Relapse interval (d)', fontsize=6); cb.ax.tick_params(labelsize=5.5)
ax = L.ax(1, 6, 6, letter='d', title='Cohort A: size of change by group')
size = dA.loc[j.index].abs().median()
for g_, c in cols.items():
    v = [size[p_] for p_ in size.index if grp[p_] == g_]
    ax.scatter([list(cols).index(g_)] * len(v), v, s=14, c=c, edgecolor='white', lw=0.3)
ax.set_xticks(range(3), ['grade 4\nIDH1-wt', 'grade 3\nIDH1-wt', 'IDH1-\nmutant'], fontsize=6); ax.set_xlim(-0.5, 2.5)
ax.set_ylabel('Median |Δ|, 699 genes')
done(L, 'Figure S2')

# ---------------------------------------------------------------- S3 all 67 genes, per patient
order = sig.sort_values('B_delta').index
L = Layout([dict(h=8.4, top=0.26, bottom=0.42, left=0.62)])
for k, (d, lab, letter) in enumerate([(dA, 'A', 'a'), (dB, 'B', 'b')]):
    ax = L.ax(0, 6 * k, 6, letter=letter, title=f'Cohort {lab}' + (' (filled, q < 0.05 in B)' if lab == 'B' else ' (no gene q < 0.05: all open)'))
    for i, g in enumerate(order):
        v = d.loc[g].values; col = UP if j.loc[g, 'B_delta'] > 0 else DN
        filled = j.loc[g, f'{lab}_q'] < 0.05
        ax.scatter(v, np.full(len(v), i) + np.linspace(-0.25, 0.25, len(v)), s=2.5, facecolor=col if filled else 'none', edgecolor=col, lw=0.3)
        ax.plot(np.median(v), i, '|', color=INK, ms=4, mew=1)
    ax.axvline(0, color=INK, lw=0.5); ax.set_ylim(-1, len(order)); ax.set_yticks(range(len(order)), order, fontsize=4.6, style='italic')
    ax.set_xlabel(f'Within-patient Δ, cohort {lab}'); ax.set_xlim(-5, 5)
done(L, 'Figure S3')

# ---------------------------------------------------------------- S4 IDH-mutant exclusion
s4 = S['S4_idh']; names = list(s4)
_dr = S['S4_dropped']; DROP = {names[0]: [], names[1]: _dr['idh_mutant'], names[2]: _dr['idh_mutant'] + _dr['grade3_wt']}
MOD4 = {'myeloid': ['CSF1R', 'SIRPA', 'MS4A6A', 'IRF8', 'DAB2', 'CMKLR1', 'AIF1', 'ITGAM', 'CD163', 'MRC1', 'TYROBP', 'FCGR3A', 'C1QA', 'C1QB'],
        'stemness': ['SOX2', 'NOTCH1', 'OLIG2', 'NES', 'POU3F2', 'DLL3', 'ASCL1'],
        'mes': ['VCAN', 'LOXL2', 'CDH11', 'FN1', 'COL1A1', 'SERPINE1', 'TGFB1', 'CD44']}
L = Layout([dict(h=2.0, top=0.26, bottom=0.62), dict(h=2.35, top=0.26, bottom=0.42, left=1.0, right=1.2)])
ax = L.ax(0, 0, 4, letter='a', title='Direction agreement')
for i, n in enumerate(names):
    ax.bar(i, s4[n]['agree_sig'], color=COH['A'], width=0.6)
    ax.text(i, s4[n]['agree_sig'] + 1.5, f"{s4[n]['agree_sig']}/67", ha='center', fontsize=6)
    ax.text(i, 4, f"p = {sci(s4[n]['binom_p'])}", ha='center', va='bottom', rotation=90, fontsize=5.8, color='white')
ax.set_xticks(range(3), ['14', '11', '10']); ax.set_ylim(0, 80); ax.set_ylabel('Same direction in A (of 67)'); ax.set_xlabel('Pairs in cohort A')
ax = L.ax(0, 4, 4, letter='b', title='Rank correlation')
for i, n in enumerate(names):
    ax.bar(i, s4[n]['rho'], color=COH['A'], width=0.6); ax.text(i, s4[n]['rho'] + 0.02, f"{s4[n]['rho']:.3f}", ha='center', fontsize=6)
ax.set_xticks(range(3), ['14', '11', '10']); ax.set_xlabel('Pairs in cohort A'); ax.set_ylim(0, 0.72); ax.set_ylabel('Spearman ρ, 699 genes')
ax = L.ax(0, 8, 4, letter='c', title='Cross-cohort prediction')
for i, n in enumerate(names):
    for off, key, tot, col in [(-0.18, 'cross_AB_pos', 17, COH['B']), (0.18, 'cross_BA_pos', s4[n]['n_pairs'], COH['A'])]:
        ax.bar(i + off, s4[n][key] / tot, 0.34, color=col)
        ax.text(i + off, s4[n][key] / tot + 0.02, f"{s4[n][key]}/{tot}", ha='center', va='bottom', fontsize=5, rotation=90)
ax.set_xticks(range(3), ['14', '11', '10']); ax.set_xlabel('Pairs in cohort A'); ax.set_ylim(0, 1.32); ax.set_ylabel('Patients with ρ > 0')
L.legend_row(0, [plt.Rectangle((0, 0), 1, 1, fc=COH['B']), plt.Rectangle((0, 0), 1, 1, fc=COH['A'])], ['A → B', 'B → A'], ncol=2, c0=8, span=4, y_off=0.0)
ax = L.ax(1, 0, 12, letter='d', title='Cohort A modules (filled, p < 0.05; open, p ≥ 0.05)')
for i, n in enumerate(names):
    for k, (m, mk) in enumerate([('myeloid', 'Myeloid'), ('stemness', 'Stemness'), ('mes', 'Mes/ECM')]):
        r = s4[n][m]; col = UP if r['median'] > 0 else DN; y = k * 4 + i
        keep_ = [c for c in dA.columns if c not in DROP[n]]; g_ = [x for x in MOD4[m] if x in dA.index]
        ax.plot(np.percentile(dA.loc[g_, keep_].mean().values, [25, 75]), [y, y], color=col, lw=0.9)
        ax.plot(r['median'], y, 'o', ms=4, mfc=col if r['p'] < 0.05 else 'white', mec=col)
        ax.text(2.35, y, f"{r['median']:+.2f}, p = {sci(r['p'])}", va='center', fontsize=5.8)
ax.set_yticks([k * 4 + i for k in range(3) for i in range(3)], [f'{m}, {p_} pairs' for m in ['Myeloid', 'Stemness', 'Mes/ECM'] for p_ in (14, 11, 10)], fontsize=5.8)
ax.invert_yaxis(); ax.axvline(0, color=INK, lw=0.5); ax.set_xlim(-2.2, 2.2)
ax.set_xlabel('Median module Δ, cohort A; bar, interquartile range of patients')
done(L, 'Figure S4')

# ---------------------------------------------------------------- S5 normalisation
s5 = pd.read_csv(f'{DER}/S5_normalisation.tsv', sep='\t')
labs = ['A\ndeposited', 'A\nupper-quartile', 'A\nquantile', 'B\nmedian-centred']
L = Layout([dict(h=2.1, top=0.26, bottom=0.50)])
for k, (col, ylab, letter, title) in enumerate([('agree', 'Same direction (%)', 'a', 'Direction agreement'), ('rho', 'Spearman ρ, 699 genes', 'b', 'Rank correlation')]):
    ax = L.ax(0, 6 * k, 6, letter=letter, title=title + ' (one cohort varied at a time)')
    v = (100 * s5.agree / s5['of']).values if col == 'agree' else s5.rho.values
    ax.bar(range(len(s5)), v, color=[COH[c] for c in s5.cohort_varied], width=0.6)
    for i, x in enumerate(v):
        ax.text(i, x + (1.5 if col == 'agree' else 0.015), f"{s5.agree[i]}/{s5['of'][i]}" if col == 'agree' else f'{x:.3f}', ha='center', fontsize=6)
    ax.set_xticks(range(len(s5)), labs, fontsize=5.8); ax.set_ylabel(ylab); ax.set_ylim(0, 112 if col == 'agree' else 0.7)
done(L, 'Figure S5')

# ---------------------------------------------------------------- S6 permutation, all nulls
nul = pd.read_csv(f'{DER}/permutation_nulls.tsv', sep='\t'); cn = pd.read_csv(f'{DER}/cross_null_allpatients.tsv', sep='\t')
pp = S['permutation']
panels = [(nul.signflip_rho, S['confirmed']['rho'], 'Sign flip: ρ', pp['signflip_rho_p'], 'ρ, 699 genes'),
          (nul.genelabel_agree, 64, 'Gene-label shuffle: agreement', pp['genelabel_agree_p'], 'same sign (of 67)'),
          (nul.genelabel_rho, S['confirmed']['rho'], 'Gene-label shuffle: ρ', pp['genelabel_rho_p'], 'ρ, 699 genes'),
          (cn.null_median_AB, S['cross']['AB_median'], 'Prediction A → B (17)', S['cross']['AB_null_p'], 'median per-patient ρ'),
          (cn.null_median_BA, S['cross']['BA_median'], 'Prediction B → A (14)', S['cross']['BA_null_p'], 'median per-patient ρ')]
L = Layout([dict(h=1.95, top=0.40, bottom=0.42), dict(h=1.95, top=0.40, bottom=0.42)])
for i, (x, obs, t, p, xl) in enumerate(panels):
    ax = L.ax(i // 3, 4 * (i % 3), 4, letter='abcde'[i], title=t,
              note=(f'observed {obs:.3f}' if obs < 1 else f'observed {obs}') + f'\np = {p:.4f}')
    ax.hist(x, bins=40, color=NEU, edgecolor='white', lw=0.2); ax.axvline(obs, color=UP, lw=1.3)
    ax.set_xlabel(xl); ax.set_ylabel('Permutations' if i % 3 == 0 else '')
done(L, 'Figure S6')

# ---------------------------------------------------------------- S7 leave-one-pair-out
lp = pd.read_csv(f'{DER}/ml_lopo_predictions.tsv', sep='\t'); s7 = pd.read_csv(f'{DER}/S7_pair_influence.tsv', sep='\t')
L = Layout([dict(h=2.0, top=0.40, bottom=0.30, left=0.40), dict(h=2.0, top=0.40, bottom=0.42)])
for i, (var, lab) in enumerate([(v, c) for v in ['sample', 'patient_centred'] for c in ['A', 'B']]):
    s_ = lp[(lp.variant == var) & (lp.cohort == lab)]
    ok = int((s_.p_recurrent > s_.p_primary).sum()); auc = S['ml'][f'{var}_LOPO_{lab}']['auc']
    ax = L.ax(0, 3 * i, 3, letter='abcd'[i], title=f"{lab}, {'sample' if var == 'sample' else 'centred'}",
              note=f'{ok}/{len(s_)} pairs ordered\nAUC {auc:.2f}')
    for _, r in s_.iterrows(): ax.plot([0, 1], [r.p_primary, r.p_recurrent], color=UP if r.p_recurrent > r.p_primary else DN, lw=0.6)
    ax.set_xticks([0, 1], ['P', 'R']); ax.set_ylim(-0.02, 1.02); ax.set_xlim(-0.2, 1.2)
    if i == 0: ax.set_ylabel('P(recurrent)')
sb = s7[s7.cohort == 'B']
for k, (col, ylab, letter, note) in enumerate([('agree', 'Agreement (fraction)', 'e', f"B-significant set: {sb['of'].min()}–{sb['of'].max()} genes"),
                                               ('rho', 'Spearman ρ, 699 genes', 'f', None)]):
    ax = L.ax(1, 6 * k, 6, letter=letter, title='One pair left out: ' + ('agreement' if col == 'agree' else 'rank correlation'), note=note, note_color=INK)
    for c in ['A', 'B']:
        s_ = s7[s7.cohort == c]; v = s_.agree / s_['of'] if col == 'agree' else s_.rho
        ax.scatter((0 if c == 'A' else 1) + np.linspace(-0.2, 0.2, len(s_)), v, s=10, c=COH[c], edgecolor='white', lw=0.3)
    ax.axhline(64 / 67 if col == 'agree' else S['confirmed']['rho'], color=UP, lw=0.8, ls=(0, (3, 2)))
    ax.set_xticks([0, 1], ['A pair left out', 'B pair left out']); ax.set_xlim(-0.5, 1.5); ax.set_ylabel(ylab)
done(L, 'Figure S7')

# ---------------------------------------------------------------- S8 resampling attenuation diagnostic
boot = pd.read_csv(f'{DER}/bootstrap.tsv', sep='\t')
DA_, DB_ = dA.loc[j.index].values, dB.loc[j.index].values
rng_ = np.random.default_rng(8); uniq, rr = [], []
for _ in range(2000):                                   # same resampling, recording distinct patients
    ia = rng_.integers(0, DA_.shape[1], DA_.shape[1]); ib = rng_.integers(0, DB_.shape[1], DB_.shape[1])
    uniq.append(len(set(ia)) + len(set(ib)))
    rr.append(stats.spearmanr(np.median(DA_[:, ia], axis=1), np.median(DB_[:, ib], axis=1))[0])
uniq, rr = np.array(uniq), np.array(rr)
L = Layout([dict(h=1.95, top=0.40, bottom=0.42), dict(h=2.1, top=0.40, bottom=0.42)])
for k, (c, obs, lab) in enumerate([('agree_sig', 64 / 67, 'Same direction, 67 genes'), ('agree_top30', 29 / 30, 'Same direction, top 30'),
                                   ('rho', S['confirmed']['rho'], 'Spearman ρ, 699 genes')]):
    ax = L.ax(0, 4 * k, 4, letter='abc'[k], title=lab, note=f'observed {obs:.3f}\nresampled median {np.median(boot[c]):.3f}')
    ax.hist(boot[c], bins=30, color=NEU, edgecolor='white', lw=0.2); ax.axvline(obs, color=UP, lw=1.3)
    ax.set_ylabel('Resamples' if k == 0 else '')
r_, p_ = stats.spearmanr(uniq, rr)
ax = L.ax(1, 0, 6, letter='d', title='Fewer distinct patients, lower ρ', note=f'all 31 patients: ρ = {S["confirmed"]["rho"]:.3f} (dashed)')
ax.scatter(uniq + rng_.uniform(-0.2, 0.2, len(uniq)), rr, s=3, facecolor='none', edgecolor=NEU, lw=0.3, rasterized=True)
med = pd.Series(rr).groupby(uniq).median(); ax.plot(med.index, med.values, color=INK, lw=1.2)
ax.axhline(S['confirmed']['rho'], color=UP, lw=1.0, ls=(0, (3, 2)))
ax.set_xlabel('Distinct patients in the resample (of 31)'); ax.set_ylabel('Resampled ρ')
ax = L.ax(1, 6, 6, letter='e', title='Distinct patients per resample',
          note=f'median {int(np.median(uniq))} of 31\nρ(distinct, resampled ρ) = {r_:.2f}', note_color=INK)
ax.hist(uniq, bins=np.arange(uniq.min() - 0.5, uniq.max() + 1.5), color=NEU, edgecolor='white', lw=0.2)
ax.set_xlabel('Distinct patients in the resample'); ax.set_ylabel('Resamples')
json.dump(dict(median_unique=float(np.median(uniq)), rho_unique_vs_rho=float(r_), p=float(p_)), open(f'{DER}/S8_resampling_diag.json', 'w'))
done(L, 'Figure S8')

# ---------------------------------------------------------------- S9 cohort-A signal outside the shared panel
gA = pd.read_csv(f'{DER}/A_genomewide.tsv', sep='\t', index_col=0); gs = pd.read_csv(f'{DER}/gsea_hallmark.tsv', sep='\t', index_col=0)
outside = ~gA.symbol.isin(j.index)
g = gs.dropna(subset=['NES_A', 'NES_B']); both = (g.q_A < 0.05) & (g.q_B < 0.05); rr9 = stats.spearmanr(g.NES_A, g.NES_B)[0]
L = Layout([dict(h=2.35, top=0.40, bottom=0.42)])
ax = L.ax(0, 0, 6, letter='a', title='Cohort A outside the shared panel',
          note=f"{int((gA.p[outside] < 0.01).sum())} of {int(outside.sum()):,} genes p < 0.01\n({round(outside.sum() * 0.01)} expected)")
ax.hist(gA.p[outside], bins=20, color=NEU, edgecolor='white', lw=0.2); ax.axhline(outside.sum() / 20, color=INK, lw=0.7, ls=(0, (3, 2)))
ax.set_xlabel('p, cohort A (genes not on the panel)'); ax.set_ylabel('Genes')
ax = L.ax(0, 6, 6, letter='b', title='Hallmark agreement',
          note=f'{len(g)} sets testable in both; ρ = {rr9:.2f}\nsame sign {int(((g.NES_A > 0) == (g.NES_B > 0)).sum())}/{len(g)}')
ax.scatter(g.NES_A[~both], g.NES_B[~both], s=12, facecolor='none', edgecolor=NEU, lw=0.6)
ax.scatter(g.NES_A[both], g.NES_B[both], s=14, c=np.where(g.NES_A[both] > 0, UP, DN), edgecolor='white', lw=0.3)
lab9 = {'E2F_TARGETS': ('e2f targets and\ng2m checkpoint', (-2.35, -1.25)), 'MTORC1_SIGNALING': ('mtorc1 signaling', (-0.7, -2.1)),
        'INTERFERON_GAMMA_RESPONSE': ('interferon γ response', (0.3, 1.95))}
for n, (t, xy) in lab9.items():
    ax.annotate(t, (g.NES_A[n], g.NES_B[n]), xytext=xy, fontsize=5.4, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.4, color=INK, shrinkA=2, shrinkB=2))
ax.set_ylim(-2.45, 2.25); ax.set_xlim(-3.15, 2.2)
ax.axhline(0, color=GRID, lw=0.6); ax.axvline(0, color=GRID, lw=0.6)
ax.set_xlabel('NES, cohort A (genome-wide)'); ax.set_ylabel('NES, cohort B (panel)')
done(L, 'Figure S9')

# ---------------------------------------------------------------- S10 module definitions
s9 = pd.read_csv(f'{DER}/S9_module_sensitivity.tsv', sep='\t')
combos = [(m, c) for m in ['Myeloid/macrophage', 'Stemness', 'Mesenchymal/ECM'] for c in ['A', 'B'] if ((s9.module == m) & (s9.cohort == c)).any()]
L = Layout([dict(h=1.95, top=0.40, bottom=0.42), dict(h=1.95, top=0.40, bottom=0.42)])
for i, (m, c) in enumerate(combos):
    s_ = s9[(s9.module == m) & (s9.cohort == c)]
    loo = s_[s_.variant.str.startswith('without')]; ext = s_[s_.variant == 'extended set'].iloc[0]
    col = UP if loo['median'].median() > 0 else DN
    ax = L.ax(i // 3, 4 * (i % 3), 4, letter='abcde'[i], title=f"{m.split('/')[0]}, cohort {c}",
              note=f"{int((loo.p < 0.05).sum())}/{len(loo)} leave-one-out p < 0.05\nextended set p = {sci(ext.p)}", note_color=INK)
    ax.scatter(np.zeros(len(loo)) + np.linspace(-0.15, 0.15, len(loo)), loo['median'], s=9, facecolor=np.where(loo.p < 0.05, col, 'white'), edgecolor=col, lw=0.6)
    ax.scatter([1], [ext['median']], s=20, marker='D', facecolor=col if ext.p < 0.05 else 'white', edgecolor=col, lw=0.8)
    ax.axhline(0, color=INK, lw=0.5); ax.set_xticks([0, 1], ['leave one\ngene out', 'extended\nset'], fontsize=6); ax.set_xlim(-0.5, 1.5)
    if i % 3 == 0: ax.set_ylabel('Median module Δ')
done(L, 'Figure S10')
print('Supplementary S1-S10 written; all letter alignment', 'OK' if OKALL else 'FAIL')
