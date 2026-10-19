#!/usr/bin/env python3
"""Figure 1 - Design and cohorts.
(a) design  (b) PCA cohort A, primary->recurrent joined per patient  (c) same, cohort B
(d) sample QC  (e) gene space  (f) cohort summary
No UMAP: 28 and 34 samples (FIG_PLAN warning)."""
import json
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.patches import FancyBboxPatch
from style import *

pcaA = pd.read_csv(f'{DER}/pca_A.tsv', sep='\t', index_col=0)
pcaB = pd.read_csv(f'{DER}/pca_B.tsv', sep='\t', index_col=0)
qcA = pd.read_csv(f'{DER}/qc_A.tsv', sep='\t', index_col=0)
qcB = pd.read_csv(f'{DER}/qc_B.tsv', sep='\t', index_col=0)
S = json.load(open(f'{DER}/summary.json'))

def partner_nn(path):
    """patients whose two samples are each other's most correlated sample (Spearman, all genes)."""
    C = pd.read_csv(path, sep='\t', index_col=0); np.fill_diagonal(C.values, -9)
    pat = {c: (c.split('_')[0] if '_' in c else (c.strip().split()[0] if ' ' in c.strip() else c.strip()[:-1])) for c in C.columns}
    best = C.idxmax(axis=1); k = 0; pats = sorted(set(pat.values()))
    for p_ in pats:
        cs = [c for c in C.columns if pat[c] == p_]
        if len(cs) == 2 and best[cs[0]] == cs[1] and best[cs[1]] == cs[0]: k += 1
    return k, len(pats)
nnA, nnB = partner_nn(f'{DER}/samplecorr_A.tsv'), partner_nn(f'{DER}/samplecorr_B.tsv')

L = Layout([dict(h=1.30, top=0.02, bottom=0.02, left=0.30),
            dict(h=2.20, top=0.40, bottom=0.62),
            dict(h=1.85, top=0.40, bottom=0.42)])
fig = L.fig

# ---------------- (a) design
ax = L.ax(0, 0, 12, letter='a', frame=False); ax.set_xlim(0, 100); ax.set_ylim(0, 27)
def box(x, y, w, h, title, lines, edge):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1', fc='white', ec=edge, lw=1.0))
    ax.text(x + 1.2, y + h - 1.0, title, fontsize=7, fontweight='bold', va='top', color=edge)
    for i, l in enumerate(lines): ax.text(x + 1.2, y + h - 4.1 - 2.75 * i, l, fontsize=6.5, va='top')
box(0.5, 14.3, 38, 12.2, 'Cohort A  GSE339484', ['Taiwan, bulk RNA-seq, 14 pairs', 'Δ = recurrent − primary,', 'per patient, inside cohort A'], COH['A'])
box(0.5, 0.6, 38, 12.2, 'Cohort B  GSE254874', ['Switzerland, NanoString IO 360, 17 pairs', 'Δ = recurrent − primary,', 'per patient, inside cohort B'], COH['B'])
box(55, 3.8, 44, 19.6, 'Compared between cohorts: Δ only', ['699 genes measured in both',
    'direction of genes significant in B', 'rank correlation of median Δ',
    'classifier trained in one cohort, tested in the other', 'module and pathway change'], INK)
for y0 in (20.4, 6.7):
    ax.annotate('', xy=(54.2, 13.6), xytext=(39.6, y0), arrowprops=dict(arrowstyle='-|>', lw=0.8, color=MUTED))
ax.text(77, 0.8, 'expression values never pooled', fontsize=6.5, color=MUTED, ha='center', style='italic')

# ---------------- (b,c) PCA with pairs joined
for k, (pc, lab, nn, key) in enumerate([(pcaA, 'A', nnA, 'pca_A'), (pcaB, 'B', nnB, 'pca_B')]):
    ax = L.ax(1, 6 * k, 6, letter='bc'[k],
              title=f'Cohort {lab}: ' + ('23,084 genes, 28 samples' if lab == 'A' else '784 genes, 34 samples'),
              note=f'own partner is the most similar sample\nin {nn[0]} of {nn[1]} patients')
    for p_, g in pc.groupby('patient'):
        q = g[g.arm == 'primary'].iloc[0]; r = g[g.arm == 'recurrent'].iloc[0]
        ax.annotate('', xy=(r.PC1, r.PC2), xytext=(q.PC1, q.PC2),
                    arrowprops=dict(arrowstyle='-|>', lw=0.6, color=NEU, shrinkA=2, shrinkB=2, mutation_scale=6))
    pr = pc[pc.arm == 'primary']; rc = pc[pc.arm == 'recurrent']
    ax.scatter(pr.PC1, pr.PC2, s=14, facecolor='white', edgecolor=COH[lab], lw=0.9, zorder=3)
    ax.scatter(rc.PC1, rc.PC2, s=14, facecolor=COH[lab], edgecolor='white', lw=0.4, zorder=4)
    ve = S[key]['var_explained']
    ax.set_xlabel(f'PC1 ({100 * ve[0]:.0f}%)'); ax.set_ylabel(f'PC2 ({100 * ve[1]:.0f}%)')
L.legend_row(1, [plt.Line2D([], [], ls='none', marker='o', ms=4, mfc='white', mec=INK, mew=0.9),
                 plt.Line2D([], [], ls='none', marker='o', ms=4, mfc=INK, mec='white'),
                 plt.Line2D([], [], color=NEU, lw=0.8)],
             ['primary', 'recurrent', 'same patient, primary → recurrent'], ncol=3, y_off=0.02)

# ---------------- (d) QC, paired (two equal sub-axes, one letter)
axs = L.subaxes(2, 0, 4, 2, letter='d', title='Sample QC, paired', gap=0.30, left=0.30)
for ax, qc, col, lab, t in [(axs[0], qcA, 'library_M', 'A', 'A: column sum (M)'),
                            (axs[1], qcB, 'median_log2cpm', 'B', 'B: median log$_2$ CPM')]:
    q = qc[qc.arm == 'primary'].set_index('patient')[col]; r = qc[qc.arm == 'recurrent'].set_index('patient')[col].loc[q.index]
    for p_ in q.index: ax.plot([0, 1], [q[p_], r[p_]], color=NEU, lw=0.6, zorder=1)
    ax.scatter(np.zeros(len(q)), q, s=9, facecolor='white', edgecolor=COH[lab], lw=0.8, zorder=3)
    ax.scatter(np.ones(len(r)), r, s=9, facecolor=COH[lab], edgecolor='white', lw=0.3, zorder=3)
    p = stats.wilcoxon(r.values - q.values).pvalue
    ax.set_xticks([0, 1], ['P', 'R']); ax.set_xlim(-0.4, 1.4)
    ax.text(0.0, 1.02, f'{t}\np = {sci(p)}', transform=ax.transAxes, ha='left', va='bottom', fontsize=6)

# ---------------- (e) gene space
ax = L.ax(2, 4, 4, letter='e', title='Gene space', note='no gene can reach q < 0.05 in A\n(min q = 0.123 at n = 14)')
bars = [('A\nexpr.', 23084, COH['A'], 0.9), ('A\nsymbol', 20373, COH['A'], 0.55), ('B\npanel', 784, COH['B'], 0.9), ('Shared', 699, INK, 0.9)]
for i, (n, v, c, al) in enumerate(bars):
    ax.bar(i, v, color=c, width=0.62, alpha=al)
    ax.text(i, v + 500, f'{v:,}' if v >= 10000 else f'{v}', ha='center', va='bottom', fontsize=6)
ax.set_xticks(range(4), [b[0] for b in bars], fontsize=6); ax.set_ylim(0, 27500); ax.set_yticks([0, 10000, 20000], ['0', '10k', '20k'])
ax.set_ylabel('Genes')

# ---------------- (f) summary table
ax = L.ax(2, 8, 4, letter='f', title='Cohort summary', frame=False, left=0.30, right=0.0)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
rows = [['', 'A', 'B'], ['Platform', 'RNA-seq', 'IO 360'], ['Pairs', '14', '17'],
        ['q < 0.05 (699)', '0', '67'], ['Relapse (d)', 'n.r.', '148–1260']]
for i, r in enumerate(rows):
    for k, t in enumerate(r):
        ax.text([0.0, 0.43, 0.72][k], 0.97 - i * 0.2, t, fontsize=6.1, va='top', transform=ax.transAxes,
                fontweight='bold' if i == 0 or k == 0 else 'normal', color=COH['AB'[k - 1]] if (i == 0 and k) else INK)
ax.plot([0, 1], [0.79, 0.79], color=MUTED, lw=0.5)
ax.text(0.0, 0.0, 'n.r., not recorded', fontsize=5.6, color=MUTED, transform=ax.transAxes)
L.check('Figure 1'); L.save('Figure1')
print('Figure 1', nnA, nnB)
