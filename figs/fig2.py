#!/usr/bin/env python3
"""Figure 2 - Primary versus recurrent change within each cohort.
(a) volcano A (genome-wide)  (b) volcano B (panel)  (c) paired slopes, top genes A
(d) paired slopes, top genes B  (e) MA plots  (f) effect-size distributions (699 genes)"""
import json, os
import numpy as np, pandas as pd
from style import *

S = json.load(open(f'{DER}/summary.json'))
gA = pd.read_csv(f'{DER}/A_genomewide.tsv', sep='\t', index_col=0)
gB = pd.read_csv(f'{DER}/B_panel784.tsv', sep='\t', index_col=0)
j = pd.read_csv(f'{DER}/TableS1_pergene_699.tsv', sep='\t', index_col=0)
H = os.path.join(HERE, '..', 'scripts')
rawA = pd.read_csv(f'{H}/GSE339484_counts_All_Samples.txt.gz', sep='\t', index_col=0)
sym = json.load(open(f'{H}/out/ens2sym.json'))
logA = np.log2(rawA + 1); logA = logA[(logA > 1).sum(axis=1) >= 14]
logA.index = [sym.get(g.split('.')[0], '') for g in logA.index]; logA = logA[logA.index != ''].groupby(level=0).mean()
b = pd.read_csv(f'{H}/GSE254874_processed.txt.gz', sep='\t', index_col=0)
def armB(c):
    c = c.strip().lower(); return c.split()[1] if ' ' in c else {'p': 'primary', 'r': 'recurrent'}[c[-1]]
def patB(c):
    c = c.strip(); return c.split()[0] if ' ' in c else c[:-1]

L = Layout([dict(h=2.15, top=0.40, bottom=0.42),
            dict(h=1.80, top=0.40, bottom=0.42),
            dict(h=1.80, top=0.40, bottom=0.42)])
ga = S['A_genomewide']

# (a) volcano A
ax = L.ax(0, 0, 6, letter='a', title='Cohort A, RNA-seq, 14 pairs',
          note=f"0 of {ga['n']:,} genes q < 0.05\n{ga['p_lt_001']} with p < 0.01 ({ga['expected']} expected)")
y = -np.log10(gA.p); hi = gA.p < 0.01
ax.scatter(gA.median_delta[~hi], y[~hi], s=2, facecolor='none', edgecolor=NEU, lw=0.3, rasterized=True)
for m, col in [(hi & (gA.median_delta > 0), UP), (hi & (gA.median_delta < 0), DN)]:
    ax.scatter(gA.median_delta[m], y[m], s=3, facecolor='none', edgecolor=col, lw=0.4, rasterized=True)
ax.axhline(-np.log10(1.22e-4), color=MUTED, lw=0.6, ls=(0, (3, 2)))
ax.text(3.95, -np.log10(1.22e-4) + 0.1, 'smallest exact p at n = 14', fontsize=5.6, color=MUTED, ha='right')
ax.set_xlim(-4, 4); ax.set_ylim(0, 4.6)
ax.set_xlabel('Median Δ (log$_2$), cohort A'); ax.set_ylabel('−log$_{10}$ p')

# (b) volcano B
sg = gB.q < 0.05
ax = L.ax(0, 6, 6, letter='b', title='Cohort B, IO 360 panel, 17 pairs',
          note=f"{int(sg.sum())} of 784 genes q < 0.05\n({int((sg & (gB.median_delta > 0)).sum())} higher, {int((sg & (gB.median_delta < 0)).sum())} lower)")
y = -np.log10(gB.p)
ax.scatter(gB.median_delta[~sg], y[~sg], s=5, facecolor='none', edgecolor=NEU, lw=0.4)
for m, col in [(sg & (gB.median_delta > 0), UP), (sg & (gB.median_delta < 0), DN)]:
    ax.scatter(gB.median_delta[m], y[m], s=8, c=col, edgecolor='white', lw=0.3)
for g, xy in {'VEGFA': (-2.25, 4.75), 'SOX2': (-0.55, 5.05), 'MS4A6A': (0.45, 5.05), 'IRF8': (1.95, 4.75)}.items():
    ax.annotate(g, (gB.loc[g, 'median_delta'], -np.log10(gB.loc[g, 'p'])), xytext=xy, fontsize=5.8, style='italic', ha='center',
                va='bottom', arrowprops=dict(arrowstyle='-', lw=0.4, color=INK, shrinkA=0, shrinkB=2))
ax.set_ylim(0, 5.5); ax.set_xlim(-2.8, 2.9)
ax.set_xlabel('Median Δ (log$_2$), cohort B'); ax.set_ylabel('−log$_{10}$ p')

# (c,d) paired slopes: three equal sub-axes per letter
pairsA = [(f'Patient{i:02d}_Primary', f'Patient{i:02d}_Recurrent') for i in range(1, 15)]
cols = list(b.columns); pats = sorted({patB(c) for c in cols})
pairsB = [([c for c in cols if patB(c) == p_ and armB(c) == 'primary'][0], [c for c in cols if patB(c) == p_ and armB(c) == 'recurrent'][0]) for p_ in pats]
topA = j.sort_values('A_p').head(3).index.tolist(); topB = j.sort_values('B_p').head(3).index.tolist()
for k, (expr, pairs, genes, letter, lab) in enumerate([(logA, pairsA, topA, 'c', 'A'), (b, pairsB, topB, 'd', 'B')]):
    axs = L.subaxes(1, 6 * k, 6, 3, letter=letter, title=f'Cohort {lab}: three smallest p, 699 shared genes', gap=0.24)
    for i, (ax, g) in enumerate(zip(axs, genes)):
        q = np.array([expr.loc[g, a] for a, _ in pairs]); r = np.array([expr.loc[g, c] for _, c in pairs]); upc = r > q
        for kk in range(len(q)): ax.plot([0, 1], [q[kk], r[kk]], color=UP if upc[kk] else DN, lw=0.55, alpha=0.8)
        ax.set_xticks([0, 1], ['P', 'R']); ax.set_xlim(-0.3, 1.3)
        ax.text(0.5, 1.03, g, transform=ax.transAxes, ha='center', va='bottom', fontsize=6.3, style='italic')
        ax.set_xlabel(f'{int(upc.sum())}/{len(q)} up', fontsize=6)
        if i == 0: ax.set_ylabel('log$_2$ expression')

# (e) MA plots, two equal sub-axes, common title
axs = L.subaxes(2, 0, 6, 2, letter='e', title='MA plots (left, cohort A; right, cohort B)', gap=0.30)
for ax, (g, xl) in zip(axs, [(gA, 'Mean log$_2$ expression'), (gB, 'Mean log$_2$ CPM')]):
    sgm = (g.q < 0.05).values
    ax.scatter(g.mean_log2[~sgm], g.median_delta[~sgm], s=1.5 if len(g) > 1000 else 4, facecolor='none', edgecolor=NEU, lw=0.3, rasterized=True)
    if sgm.any():
        ax.scatter(g.mean_log2[sgm], g.median_delta[sgm], s=6, c=np.where(g.median_delta[sgm] > 0, UP, DN), edgecolor='white', lw=0.2)
    ax.axhline(0, color=INK, lw=0.5); ax.set_xlabel(xl); ax.set_ylim(-3.5, 3.5)
axs[0].set_ylabel('Median Δ (log$_2$)')

# (f) effect-size distribution, each cohort on its own scale
ax = L.ax(2, 6, 6, letter='f', title='Median Δ per gene, own scale (not pooled)')
bins = np.linspace(-2.5, 2.5, 51)
x0, y0, w0, h0 = L.box(2, 6, 6)
for k, (lab, col) in enumerate([('A', COH['A']), ('B', COH['B'])]):
    v = j[f'{lab}_delta']; ax.hist(v, bins=bins, histtype='step', color=col, lw=1.1)
    L.fig.text(x0 / L.W, (y0 + h0 + 0.40 - 0.02 - LINE * (k + 1)) / L.H, f'cohort {lab}: SD {v.std():.2f}, {100 * (v > 0).mean():.0f}% of genes > 0',
               fontsize=6.3, color=col, fontweight='bold', va='top')
ax.set_xlabel('Median Δ (log$_2$), 699 shared genes'); ax.set_ylabel('Genes')
L.check('Figure 2'); L.save('Figure2')
print('Figure 2 written', topA, topB)
