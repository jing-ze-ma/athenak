#!/usr/bin/env python3
"""cost_pie.py: draw cost_pie.png / .svg from cost_modules.csv (same directory).  Style: mesa2.mplstyle (idort_style.py)."""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, Wedge
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('IDORT_STYLE', os.path.join(HERE, '..', '..', '..', 'tests', 'vet_clean', 'notebooks')))
import idort_style
idort_style.use()
COL = {'rad-vet': '#0072B2', 'rays': '#E69F00', 'hydro': '#009E73', 'others': '#CC79A7'}   # Okabe & Ito
ORDER = ['rad-vet', 'rays', 'hydro', 'others']                                             # counter-clockwise from the top
TITLE = {'bsg': 'M20 BSG', 'm35mms': 'M35 MMS', 'm35zams': 'M35 ZAMS', 't9l62': 'M56 LBV (radiative IC) *'}
light = lambda c, f=0.55: tuple(1 - (1 - x)*f for x in to_rgb(c))
D = {}
for r in csv.DictReader(open(os.path.join(HERE, 'cost_modules.csv'))):
    d = D.setdefault(r['star'], {'cells': int(r['cells']), 'p': {}}); d['p'][r['module'], r['part']] = float(r['percent'])
pct = lambda v: f'{v:.0f}%' if v >= 0.95 else f'{v:.1f}%'
at = lambda o, r, a: (o[0] + r*np.cos(np.radians(a)), o[1] + r*np.sin(np.radians(a)))
R0, R1, R2, A0, A1 = 0.42, 0.79, 0.99, 1.09, 1.16
fig, axs = plt.subplots(2, 2, figsize=(14, 14.4))
for ax, (st, d) in zip(axs.flat, D.items()):
    p = d['p']; a = 90.0
    for m in ORDER:
        v = sum(x for (mm, _), x in p.items() if mm == m); w = 3.6*v; mid = a + w/2; o = at((0, 0), 0.05, mid)
        ax.add_patch(Wedge(o, R1, a, a + w, width=R1 - R0, fc=light(COL[m]), ec='none'))
        ax.add_patch(Wedge(o, R2, a, a + w, width=R2 - R1, fc=COL[m], ec='none'))
        if v > 2:
            ax.text(*at(o, (R1 + R2)/2, mid), pct(v), ha='center', va='center', fontsize=15 if v > 10 else 12,
                    color='white' if m == 'rad-vet' else 'black')
            nm = {'others': 'others'}.get(m, m)
            ax.text(*at(o, (R0 + R1)/2 + (0.02 if v < 10 else 0), mid), nm, ha='center', va='center', fontsize=(14.5 if m == 'rad-vet' else 16) if v > 10 else 13,
                    weight='bold' if m == 'rad-vet' else 'normal', rotation=0 if v > 10 else (mid - 180 if 90 < mid < 270 else mid if mid < 90 else mid - 360))
        else: ax.text(*at(o, A1 + 0.05, mid), f'{m} {pct(v)}', ha='right', va='bottom', fontsize=13)
        parts = [('calculation', dict(fc=COL[m], ec='none')), ('communication', dict(fc=light(COL[m], 0.3), ec=COL[m], hatch='/////', lw=0.8)),
                 ('rest', dict(fc='white', ec=COL[m], lw=1.2))]
        b = a
        for part, kw in parts:
            x = p.get((m, part), 0) if part != 'rest' else p.get((m, 'unattributed'), 0) + p.get((m, 'GPU idle'), 0)
            if x < 0.05: continue
            s = 3.6*x; g = min(0.6, s/6)
            ax.add_patch(Wedge(o, A1, b + g, b + s - g, width=A1 - A0, **kw))
            if x >= 0.95 and v > 2:
                q = at(o, A1 + 0.12, b + s/2)
                ax.text(*q, pct(x), ha='center', va='center', fontsize=13)
            b += s
        a += w
    ax.text(0, 0, f'{d["cells"]/1e6:.1f} M cells\n4 H200\n(1 node)', ha='center', va='center', fontsize=17, linespacing=1.35)
    ax.set_title(TITLE[st], fontsize=22, pad=8); ax.set_xlim(-1.42, 1.42); ax.set_ylim(-1.34, 1.34); ax.set_aspect('equal'); ax.axis('off')
H = [Patch(fc='0.35', ec='none', label='calculation'), Patch(fc='0.85', ec='0.35', hatch='/////', label='communication'),
     Patch(fc='white', ec='0.35', lw=1.2, label='GPU idle without MPI, unattributed')]
fig.legend(handles=H, loc='lower center', ncol=3, fontsize=17, bbox_to_anchor=(0.5, 0.012), handlelength=2.6)
fig.text(0.5, 0.0, '* fenced kernel timer (another measurement than the other three; see the text)', ha='center', fontsize=13)
fig.subplots_adjust(left=0.01, right=0.99, top=0.965, bottom=0.055, wspace=0.0, hspace=0.06)
fig.savefig(os.path.join(HERE, 'cost_pie.svg')); fig.savefig(os.path.join(HERE, 'cost_pie.png'), dpi=110)
