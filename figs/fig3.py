#!/usr/bin/env python3
"""Figure 3 - Reproducibility across cohorts and platforms (the paper's main figure).
(a) 64/67 same direction  (b) top 30, 29/30  (c) 699-gene scatter, rho in panel
(d) CAT curves  (e) classifier trained in one cohort, tested in the other (ROC, AUC)
(f) sign-flip permutation null, agreement  (g) sign-flip permutation null, rho"""
import json
import numpy as np, pandas as pd
from scipy import stats
from style import *

S = json.load(open(f'{DER}/summary.json'))
j = pd.read_csv(f'{DER}/TableS1_pergene_699.tsv', sep='\t', index_col=0)
cat = pd.read_csv(f'{DER}/cat.tsv', sep='\t')
roc = pd.read_csv(f'{DER}/ml_roc.tsv', sep='\t')
nul = pd.read_csv(f'{DER}/permutation_nulls.tsv', sep='\t')
sig = j[j.B_q < 0.05]; top = j.sort_values('B_q').head(30)
pp = S['permutation']; ml = S['ml']

L = Layout([dict(h=2.05, top=0.40, bottom=0.62),
            dict(h=2.15, top=0.40, bottom=0.62),
            dict(h=1.90, top=0.40, bottom=0.42)])

def dirbars(ax, sub):
    up = sub[sub.B_delta > 0]; dn = sub[sub.B_delta < 0]
    for i, (s_, col) in enumerate([(up, UP), (dn, DN)]):
        ok = int(((s_.A_delta > 0) == (s_.B_delta > 0)).sum()); bad = len(s_) - ok
        ax.bar(i, ok, color=col, width=0.6)
        ax.bar(i, bad, bottom=ok, facecolor='white', edgecolor=col, lw=0.9, width=0.6)
        ax.text(i, len(s_) + max(len(up), len(dn)) * 0.03, f'{ok}/{len(s_)}', ha='center', va='bottom', fontsize=6.3)
    ax.set_xticks([0, 1], ['higher\nin B', 'lower\nin B']); ax.set_ylabel('Genes')
    ax.set_ylim(0, max(len(up), len(dn)) * 1.18); ax.set_xlim(-0.6, 1.6)

ax = L.ax(0, 0, 4, letter='a', title='67 genes with q < 0.05 in B',
          note=f"64/67 = 95.5% same sign in A\nbinomial p = {sci(S['confirmed']['binom_p'])}"); dirbars(ax, sig)
ax = L.ax(0, 4, 4, letter='b', title='Top 30 in B',
          note=f"29/30 = 96.7% same sign in A\nbinomial p = {sci(S['confirmed']['binom_p30'])}"); dirbars(ax, top)
L.legend_row(0, [plt.Rectangle((0, 0), 1, 1, fc=INK), plt.Rectangle((0, 0), 1, 1, fc='white', ec=INK)],
             ['same sign in A', 'opposite sign in A'], ncol=2, c0=0, span=8, y_off=0.0)

# (c) scatter
ax = L.ax(0, 8, 4, letter='c', title='699 shared genes',
          note=f"ρ = {S['confirmed']['rho']:.3f}\np = {sci(S['confirmed']['rho_p'])}")
ns = j[j.B_q >= 0.05]
ax.scatter(ns.A_delta, ns.B_delta, s=5, facecolor='none', edgecolor=NEU, lw=0.5, zorder=1)
for s_, col in [(sig[sig.B_delta > 0], UP), (sig[sig.B_delta < 0], DN)]:
    ax.scatter(s_.A_delta, s_.B_delta, s=9, c=col, edgecolor='white', lw=0.3, zorder=3)
b1, b0 = np.polyfit(j.A_delta, j.B_delta, 1); xx = np.linspace(j.A_delta.min(), j.A_delta.max(), 10)
ax.plot(xx, b0 + b1 * xx, color=INK, lw=0.8, zorder=2)
for g, xy in {'ANGPT1': (-0.85, 2.3), 'ANLN': (1.75, -1.75), 'VEGFA': (0.95, -2.45)}.items():
    ax.annotate(g, (j.loc[g, 'A_delta'], j.loc[g, 'B_delta']), xytext=xy, fontsize=5.8, style='italic', ha='center', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.4, color=INK, shrinkA=3, shrinkB=2))
ax.axhline(0, color=GRID, lw=0.6, zorder=0); ax.axvline(0, color=GRID, lw=0.6, zorder=0)
ax.set_ylim(-2.75, 2.75); ax.set_xlim(-1.8, 2.3)
ax.set_xlabel('Median Δ, cohort A'); ax.set_ylabel('Median Δ, cohort B')

# (d) CAT
c30, c100 = cat[cat.k == 30].iloc[0], cat[cat.k == 100].iloc[0]
ax = L.ax(1, 0, 6, letter='d', title='Concordance at the top',
          note=f'top 30: {c30.up:.2f} higher, {c30.down:.2f} lower (chance {c30.expected:.2f})\n'
               f'top 100: {c100.up:.2f} higher, {c100.down:.2f} lower (chance {c100.expected:.2f})', note_color=INK)
ln1, = ax.plot(cat.k, cat.up, color=UP, lw=1.2); ln2, = ax.plot(cat.k, cat.down, color=DN, lw=1.2)
ln3, = ax.plot(cat.k, cat.expected, color=MUTED, lw=0.8, ls=(0, (3, 2)))
ax.set_xlabel('Top k genes in each cohort (signed −log$_{10}$ p)'); ax.set_ylabel('Shared fraction of top k')
ax.set_ylim(0, 1); ax.set_xlim(1, 300)
L.legend_row(1, [ln1, ln2, ln3], ['higher at recurrence', 'lower at recurrence', 'chance'], ncol=3, c0=0, span=6, y_off=0.0)

# (e) ROC
ax = L.ax(1, 6, 6, letter='e', title='Primary vs recurrent, elastic net',
          note='solid: trained in one cohort, tested in the other\ndotted: leave-one-pair-out within a cohort', note_color=INK)
v = 'sample'
for curve, col, ls in [('A->B', COH['B'], '-'), ('B->A', COH['A'], '-'), ('LOPO A', COH['A'], ':'), ('LOPO B', COH['B'], ':')]:
    r = roc[(roc.variant == v) & (roc.curve == curve)]
    ax.plot(r.fpr, r.tpr, color=col, lw=1.1 if ls == '-' else 0.9, ls=ls)
ax.plot([0, 1], [0, 1], color=NEU, lw=0.6, ls=(0, (3, 2)))
ax.set_xlabel('False-positive rate'); ax.set_ylabel('True-positive rate'); ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
for i, (curve, col, ls) in enumerate([('A->B', COH['B'], '—'), ('B->A', COH['A'], '—'), ('LOPO_A', COH['A'], '···'), ('LOPO_B', COH['B'], '···')]):
    name = {'A->B': 'train A, test B', 'B->A': 'train B, test A', 'LOPO_A': 'within A', 'LOPO_B': 'within B'}[curve]
    ax.text(0.98, 0.05 + 0.1 * (3 - i), f'{ls} {name}: AUC {ml[f"{v}_{curve}"]["auc"]:.2f}', transform=ax.transAxes, ha='right',
            fontsize=6, color=col, fontweight='bold' if ls == '—' else 'normal')

# (f) sign-flip null, agreement
ax = L.ax(2, 0, 6, letter='f', title='Null: random sign flip of each cohort-A patient',
          note=f'observed 64/67, permutation p = {pp["signflip_agree_p"]:.3f}\n(sign flips keep gene–gene correlation)')
x = nul.signflip_agree.values
ax.hist(x, bins=np.arange(x.min() - 0.5, 68), color=NEU, edgecolor='white', lw=0.3)
ax.axvline(64, color=UP, lw=1.4)
ax.set_xlabel('Genes with the same sign in A (of 67)'); ax.set_ylabel('Permutations')

# (g) sign-flip null, rho
obs = S['confirmed']['rho']; z_, se_ = np.arctanh(obs), 1 / np.sqrt(len(j) - 3)
lo_, hi_ = np.tanh(z_ - 1.96 * se_), np.tanh(z_ + 1.96 * se_)
ax = L.ax(2, 6, 6, letter='g', title='Null: random sign flip, rank correlation',
          note=f'observed ρ = {obs:.3f}, permutation p = {pp["signflip_rho_p"]:.3f}\nFisher z 95% CI {lo_:.2f}–{hi_:.2f} (genes independent)')
x = nul.signflip_rho.values
ax.hist(x, bins=50, color=NEU, edgecolor='white', lw=0.3)
top_ = ax.get_ylim()[1]; ax.set_ylim(0, top_ * 1.15)
ax.annotate('', xy=(obs, top_ * 0.03), xytext=(obs, top_ * 0.75), arrowprops=dict(arrowstyle='-|>', lw=1.3, color=UP, mutation_scale=9))
ax.set_xlabel('Spearman ρ of median Δ, 699 genes'); ax.set_ylabel('Permutations')
L.check('Figure 3'); L.save('Figure3')
print('Figure 3 written')
