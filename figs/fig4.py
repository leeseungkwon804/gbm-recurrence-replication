#!/usr/bin/env python3
"""Figure 4 - Modules and pathways.
(a) module forest, both cohorts  (b,c) module genes x patients, within-patient change
(d) Hallmark preranked GSEA  (e,f) paired module scores, primary vs recurrent"""
import json, os
import numpy as np, pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from style import *

mod = json.load(open(f'{DER}/modules.json'))
ms = pd.read_csv(f'{DER}/module_sample_scores.tsv', sep='\t')
gs = pd.read_csv(f'{DER}/gsea_hallmark.tsv', sep='\t', index_col=0)
H = os.path.join(HERE, '..', 'scripts')
dA = pd.read_csv(f'{H}/out/A_deltas.tsv', sep='\t', index_col=0)
sym = json.load(open(f'{H}/out/ens2sym.json'))
dA.index = [sym.get(g.split('.')[0], '') for g in dA.index]; dA = dA[dA.index != ''].groupby(level=0).mean()
dB = pd.read_csv(f'{H}/out/B_deltas.tsv', sep='\t', index_col=0)
DIV = LinearSegmentedColormap.from_list('div', [DN, '#f2f3f5', UP])
MODS = ['Myeloid/macrophage', 'Stemness', 'Mesenchymal/ECM']

L = Layout([dict(h=1.75, top=0.30, bottom=0.40),
            dict(h=3.05, top=0.26, bottom=0.30),
            dict(h=2.05, top=0.40, bottom=0.42)])

# (a) forest; values and the not-tested note sit in the right part of the axes, clear of the data
ax = L.ax(0, 0, 12, letter='a', title='Module change per cohort (filled, p < 0.05; open, p ≥ 0.05)', left=1.05)
yl = []; i = 0
for m in MODS:
    for lab in ['A', 'B']:
        r = mod.get(f'{m}|{lab}')
        if r is None:
            ax.text(2.6, i, 'not tested (2 genes on panel)', fontsize=6, color=MUTED, va='center')
        else:
            sigf = r['p'] < 0.05; col = UP if r['median'] > 0 else DN
            ax.plot(np.percentile(r['values'], [25, 75]), [i, i], color=col, lw=1.0)
            ax.plot(r['median'], i, 'o', ms=5, mfc=col if sigf else 'white', mec=col, mew=1.1)
            ax.text(2.6, i, f"{r['median']:+.2f}   p = {sci(r['p'])}   {r['pct_up']:.0f}% of patients higher", va='center', fontsize=6.2,
                    color=col if sigf else INK, fontweight='bold' if sigf else 'normal')
        yl.append(f'{m}, {lab}'); i += 1
ax.axvline(0, color=INK, lw=0.5); ax.set_yticks(range(len(yl)), yl); ax.invert_yaxis(); ax.set_xlim(-2.2, 6.2)
ax.set_xticks([-2, -1, 0, 1, 2]); ax.spines['bottom'].set_bounds(-2.2, 2.2)
ax.set_xlabel('Median module Δ (log$_2$); bar, interquartile range of patients', loc='left')

# (b,c) heatmaps with equal cell size; cohort B has fewer module genes on the panel, so its block is shorter
def rows_for(lab):
    rows, seps = [], []
    for m in MODS:
        g = mod.get(f'{m}|{lab}')
        if g is None: continue
        rows += g['genes']; seps.append(len(rows))
    return rows, seps
rA, sA = rows_for('A'); rB, sB = rows_for('B')
x0, y0, w0, h0 = L.box(1, 0, 6, left=0.62)
cell_h = h0 / len(rA)
for k, (d, rows, seps, lab, letter) in enumerate([(dA, rA, sA, 'A', 'b'), (dB, rB, sB, 'B', 'c')]):
    x, y, w, h = L.box(1, 6 * k, 6, left=0.62)
    hh = cell_h * len(rows)
    ax = L.fig.add_axes([x / L.W, (y + h - hh) / L.H, w / L.W, hh / L.H])
    im = ax.imshow(d.loc[rows].values, aspect='auto', cmap=DIV, norm=TwoSlopeNorm(0, -3, 3), interpolation='nearest')
    for s_ in seps[:-1]: ax.axhline(s_ - 0.5, color='white', lw=1.5)
    ax.set_yticks(range(len(rows)), rows, fontsize=5.4, style='italic'); ax.set_xticks([])
    ax.set_xlabel(f'Patients, cohort {lab} (n = {d.shape[1]})', fontsize=6.5)
    L.letter(1, 6 * k, 6, letter); L.head(1, 6 * k, 6, f'Within-patient Δ, cohort {lab}', left=0.62)
# one colour bar for both heatmaps, in the free space under (c)
x, y, w, h = L.box(1, 6, 6, left=0.62)
cax = L.fig.add_axes([(x + 0.25 * w) / L.W, (y + 0.35) / L.H, 0.5 * w / L.W, 0.08 / L.H])
cb = L.fig.colorbar(im, cax=cax, orientation='horizontal'); cb.set_label('Within-patient Δ (log$_2$), both heatmaps', fontsize=6)
cb.ax.tick_params(labelsize=5.8)

# (d) GSEA
ax = L.ax(2, 0, 6, letter='d', title='Hallmark preranked GSEA', left=1.15, head_left=0.30, note_bold=False,
          note='bars: cohort A, genome-wide; diamonds: cohort B, panel\nfilled, q < 0.05; open, q ≥ 0.05', note_color=INK)
g = gs.dropna(subset=['NES_A']).sort_values('q_A').head(10).sort_values('NES_A')
for k, (name, r) in enumerate(g.iterrows()):
    col = UP if r.NES_A > 0 else DN; sigf = r.q_A < 0.05
    ax.barh(k, r.NES_A, color=col if sigf else 'white', edgecolor=col, lw=0.8, height=0.6)
    if not np.isnan(r.get('NES_B', np.nan)):
        ax.plot(r.NES_B, k, 'D', ms=3.8, mfc=COH['B'] if r.q_B < 0.05 else 'white', mec=COH['B'], mew=0.9)
ax.set_yticks(range(len(g)), [n.replace('_', ' ').lower() for n in g.index], fontsize=5.8); ax.axvline(0, color=INK, lw=0.5)
ax.set_xlabel('Normalised enrichment score')

# (e,f) paired module scores
for k, (lab, letter, c0) in enumerate([('A', 'e', 6), ('B', 'f', 9)]):
    sub = ms[ms.cohort == lab]; mods = [m for m in MODS if m in set(sub.module)]
    short = {'Myeloid/macrophage': 'Mye', 'Stemness': 'Stem', 'Mesenchymal/ECM': 'ECM'}
    counts = []
    for m in mods:
        s_ = sub[sub.module == m]; counts.append(f'{short[m]} {int((s_.recurrent > s_.primary).sum())}/{len(s_)}')
    ax = L.ax(2, c0, 3, letter=letter, title=f'Cohort {lab}, P → R', note='under each module:\npatients higher at R', note_color=INK, note_bold=False,
              left=0.42 if lab == 'A' else 0.30)
    for i, m in enumerate(mods):
        s_ = sub[sub.module == m]; z = lambda v: v - s_[['primary', 'recurrent']].values.mean()
        for _, r in s_.iterrows():
            ax.plot([3 * i, 3 * i + 1], [z(r.primary), z(r.recurrent)], color=UP if r.recurrent > r.primary else DN, lw=0.45, alpha=0.75)
        ax.boxplot([z(s_.primary), z(s_.recurrent)], positions=[3 * i, 3 * i + 1], widths=0.6, showfliers=False,
                   medianprops=dict(color=INK, lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.set_xticks([3 * i + 0.5 for i in range(len(mods))], [f'{short[m]}\n{c.split()[1]}' for m, c in zip(mods, counts)], fontsize=6)
    ax.set_xlim(-0.8, 3 * len(mods) - 1.2); ax.set_ylim(-2.4, 2.4)
    if lab == 'A': ax.set_ylabel('Module score, centred')
L.check('Figure 4'); L.save('Figure4')
print('Figure 4 written')
