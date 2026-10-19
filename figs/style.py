"""Shared figure style and layout (HQ_c5_FIG_FIX_20261014).

Style after the Jung IJMS figures: Arial; bold lower-case bracketed panel letters; colour carries
meaning (red higher at recurrence, blue lower, grey context); non-significant points hollow; key
numbers written in the figure.

Layout rules:
  * the figure is 180 mm wide (or narrower if the 230 mm height limit binds), divided into 12 equal columns; a panel cell spans
    an integer number of columns;
  * every row has one cell height; inside a cell the left / top / bottom paddings are the same for
    all cells of the row, so axes in a row share their top and bottom, and equal spans give equal
    axes widths;
  * panel letters are placed in figure coordinates at the top-left corner of the cell, so letters in
    a row share y and letters starting in the same column share x;
  * titles sit left-aligned above the axes and numeric notes on the line(s) below the title, outside
    the data area.
"""
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager as fm
fm._load_fontmanager(try_read_cache=False)
import matplotlib.pyplot as plt

W0 = 7857 / 1440                        # 5.456 in, width the row heights were designed at
W = 180 / 25.4                          # 7.087 in = 180 mm, inserted at native size (HQ c5fmt 20261015)
HMAX = 230 / 25.4                       # figure height limit, legend must fit on the same page
UP, DN, NEU, INK, GRID, MUTED = '#b2182b', '#2166ac', '#b8bcc4', '#1a1d22', '#e9ebef', '#6b7280'
COH = {'A': '#1b7837', 'B': '#8a5a9e'}  # cohort identity (validated pair)
HERE = os.path.dirname(os.path.abspath(__file__))
DER = os.path.join(HERE, '..', 'data', 'derived')

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'axes.linewidth': 0.6, 'axes.edgecolor': MUTED, 'xtick.major.width': 0.6,
    'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.color': INK,
    'ytick.color': INK, 'text.color': INK, 'axes.labelcolor': INK, 'axes.spines.top': False,
    'axes.spines.right': False, 'legend.frameon': False, 'savefig.dpi': 600,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'axes.titlepad': 0})

LINE = 0.115                            # one 7-pt text line, inches


def sci(p):
    """p-value as text, scientific below 1e-3."""
    if p >= 0.1:
        return f'{p:.2f}'
    if p >= 0.001:
        return f'{p:.3f}'
    m, e = f'{p:.1e}'.split('e')
    return f'{m} × 10$^{{{int(e)}}}$'


class Layout:
    """rows: list of dicts with h (cell height, in) and optional top/bottom/left paddings (in)."""

    def __init__(self, rows, ncol=12, margin=0.04, gutter=0.12, vgap=0.06):
        # row heights scale with the canvas width; paddings (text room, in pt-sized inches) stay as designed.
        # If the height limit binds, the scale is reduced so the figure is exactly HMAX tall and the width follows.
        fixed = vgap * (len(rows) - 1) + 2 * margin; base = sum(r['h'] for r in rows)
        k = W / W0
        if k * base + fixed > HMAX: k = (HMAX - 0.002 - fixed) / base
        self.k, self.W = k, min(W, W0 * k)
        rows = [dict(r, h=r['h'] * k) for r in rows]
        self.rows, self.ncol, self.m, self.g, self.vgap = rows, ncol, margin, gutter, vgap
        self.H = sum(r['h'] for r in rows) + fixed
        self.fig = plt.figure(figsize=(self.W, self.H))
        self.colw = (self.W - 2 * margin - (ncol - 1) * gutter) / ncol
        self.letters = []

    def cell(self, r, c0, span):
        x = self.m + c0 * (self.colw + self.g)
        w = span * self.colw + (span - 1) * self.g
        top = self.H - self.m - sum(self.rows[i]['h'] + self.vgap for i in range(r))
        return x, top - self.rows[r]['h'], w, self.rows[r]['h']

    def pads(self, r, left=None, right=None, top=None, bottom=None):
        R = self.rows[r]
        return (R.get('left', 0.45) if left is None else left, R.get('right', 0.06) if right is None else right,
                R.get('top', 0.30) if top is None else top, R.get('bottom', 0.42) if bottom is None else bottom)

    def box(self, r, c0, span, **kw):
        x, y, w, h = self.cell(r, c0, span)
        L, Rt, T, B = self.pads(r, **kw)
        return x + L, y + B, w - L - Rt, h - T - B

    def _rect(self, x, y, w, h):
        return [x / self.W, y / self.H, w / self.W, h / self.H]

    def letter(self, r, c0, span, s):
        x, y, w, h = self.cell(r, c0, span)
        fx, fy = x / self.W, (y + h) / self.H
        self.fig.text(fx, fy, f'({s})', fontsize=8.5, fontweight='bold', ha='left', va='top')
        self.letters.append((s, r, c0, round(fx, 4), round(fy, 4)))

    def head(self, r, c0, span, title=None, note=None, note_color=UP, left=None, note_bold=True):
        """title on the first line above the axes, note lines below it (still above the axes)."""
        x, y, w, h = self.box(r, c0, span, left=left)
        top = y + h
        T = self.pads(r)[2]
        if title:
            self.fig.text(x / self.W, (top + T - 0.02) / self.H, title, fontsize=7, ha='left', va='top')
        if note:
            self.fig.text(x / self.W, (top + T - 0.02 - (LINE if title else 0)) / self.H, note, fontsize=6.3,
                          ha='left', va='top', color=note_color, fontweight='bold' if note_bold else 'normal',
                          linespacing=1.15)

    def ax(self, r, c0, span, letter=None, title=None, note=None, note_color=UP, frame=True, head_left=None, note_bold=True, **pad):
        x, y, w, h = self.box(r, c0, span, **pad)
        a = self.fig.add_axes(self._rect(x, y, w, h))
        if not frame: a.axis('off')
        if letter: self.letter(r, c0, span, letter)
        self.head(r, c0, span, title, note, note_color, left=pad.get('left') if head_left is None else head_left, note_bold=note_bold)
        return a

    def subaxes(self, r, c0, span, n, letter=None, title=None, note=None, gap=0.08, sharey=False, **pad):
        """n equal axes filling the axes box of one cell (one letter)."""
        x, y, w, h = self.box(r, c0, span, **pad)
        sw = (w - gap * (n - 1)) / n; out = []
        for i in range(n):
            a = self.fig.add_axes(self._rect(x + i * (sw + gap), y, sw, h), sharey=out[0] if (sharey and out) else None)
            out.append(a)
        if letter: self.letter(r, c0, span, letter)
        self.head(r, c0, span, title, note, left=pad.get('left'))
        return out

    def legend_row(self, r, handles, labels, ncol, c0=0, span=12, y_off=0.04, **kw):
        """legend in the bottom padding of row r, centred on the given span."""
        x, y, w, h = self.cell(r, c0, span)
        self.fig.legend(handles, labels, loc='lower center', ncol=ncol, bbox_to_anchor=((x + w / 2) / self.W, (y + y_off) / self.H),
                        handlelength=1.6, columnspacing=1.2, handletextpad=0.4, **kw)

    def check(self, name):
        """letters in a row share y; letters starting in the same column share x."""
        ok = True
        print(f'{name}: {self.W * 25.4:.1f} x {self.H * 25.4:.1f} mm, row scale {self.k:.3f}; panel letters (figure coordinates)')
        for s, r, c0, fx, fy in self.letters:
            print(f'   ({s}) row {r} col {c0:2d}  x = {fx:.4f}  y = {fy:.4f}')
        for key, idx in (('row', 1), ('col', 2)):
            groups = {}
            for L in self.letters: groups.setdefault(L[idx], set()).add(L[4] if key == 'row' else L[3])
            for g, vals in groups.items():
                if len(vals) > 1: ok = False; print(f'   MISALIGNED {key} {g}: {vals}')
        print(f'   alignment: {"OK" if ok else "FAIL"}')
        return ok

    def save(self, stem, folder=HERE):
        for ext in ('pdf', 'png'):
            self.fig.savefig(os.path.join(folder, f'{stem}.{ext}'), dpi=600, facecolor='white')
        plt.close(self.fig)
