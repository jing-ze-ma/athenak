#!/usr/bin/env python3
"""Overplot published 3-D RHD massive-star ENVELOPE simulations on the HR diagram
of Marchant & Bodensteiner 2024 (ARA&A 62, 21), Figure 2.

Background: figures/HR.pdf from the arXiv source of arXiv:2311.01865, rasterised
with ghostscript at 300 dpi, cropped to the axes frame, then (a) the f_bin legend text
block masked out with white over full-res-crop pixels x=[8,975], y=[535,1570] of
HR_fig2_axes.png and (b) desaturated to grey, out = 255 - (255 - luminance)*0.80,
giving HR_fig2_axes_grey.png. No labels are drawn on the points; identities are in
sims_3d.csv / the page table.

Pixel -> data calibration (major ticks, full-resolution 6250x4063 raster):
  x (log10 Teff): ticks 5.0,4.5,4.0,3.5 at px 1710.5, 3080.5, 4450.5, 5820.5
      log10 Teff = 5.624270070 - 3.649635e-4 * px      (residual < 1e-15 dex)
  y (log10 L/Lsun): ticks 6,5,4,3 at px 1056.5, 1778.5, 2500.5, 3223.0
      log10 L    = 7.462853890 - 1.3847538e-3 * py     (max residual 2.8e-4 dex)
Cropped raster HR_fig2_axes.png covers full-res pixels x=[343.5,6145.5],
y=[410.5,3652.5] -> the EXTENT constants below.
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))

AX, BX = 5.624270070, -3.649635036e-4
AY, BY = 7.462853890, -1.38475382e-3
X_LEFT, X_RIGHT = AX + BX*343.5, AX + BX*6145.5
Y_TOP, Y_BOT = AY + BY*410.5, AY + BY*3652.5

# our own added region: the stripped-He / sdO strip below the classical WR box.
# Extent chosen from the Goetberg et al. 2018 stripped-star grid (2-7 Msun,
# log L 3.0-4.6, Teff 50-100 kK) and the Woosley 2019 He stars up to ~10 Msun
# (log L ~5.0, Teff ~125 kK); MB24 has no such box (cWR covers log L > 5.5).
# A parallelogram tilted along the helium main sequence, like the tilted O/B/OBe
# bands of MB24 follow the hydrogen MS.  The hot edge follows the Woosley (2019)
# He-star ZAMS (3 Msun: log L 3.91, log Teff 4.91; 5 Msun: 4.47, 4.98;
# 8 Msun: 4.90, 5.04 -> ~7 dex in L per dex in Teff); the cool edge is offset
# 0.30 dex cooler, for stripped stars that keep a residual H layer
# (Goetberg et al. 2018).  (log Teff, log L) corners:
HE_POLY = [(4.75, 2.8), (5.11, 5.3), (4.81, 5.3), (4.45, 2.8)]
HE_MID = ((4.60, 2.8), (4.96, 5.3))     # for the label angle
HE_LABEL_AT = (4.80, 4.10)

# on black, the blue and the purple need more luminous variants
DARK_COLOR = {'MPI-AMRVAC': '#5b9bff', 'AREPO': '#b57bff',
              'AthenaK': '#ff7a5c'}

BG_FADE = 0.5         # background brightness in the highlight versions
EDGE = '#33302a'       # thin dark marker edge (light version)
EDGE_DARK = '#f2efe8'  # thin light marker edge (dark version)

CODE_STYLE = {
    'Athena++':    dict(color='#f2c400', marker='o'),
    'AthenaK':     dict(color='#e8452c', marker='*'),
    'MPI-AMRVAC':  dict(color='#1f5fd6', marker='s'),
    'AREPO':       dict(color='#7d3ac1', marker='D'),
    'other':       dict(color='#555555', marker='v'),
}


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if not r['tag'] or r['tag'].startswith('#'):
                continue
            try:
                r['x'] = np.log10(float(r['Teff_K']))
                r['y'] = float(r['logL_Lsun'])
            except (ValueError, TypeError):
                r['x'] = r['y'] = None
            rows.append(r)
    return rows


def draw(rows, background, out, dark=False, transparent=False, codes=None,
         highlight=None):
    fig, ax = plt.subplots(figsize=(11.5, 7.6))
    fg = 'white' if dark else 'black'
    if dark:
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')
        for sp in ax.spines.values():
            sp.set_color(fg)
        ax.tick_params(colors=fg, which='both')
    if background:
        img = plt.imread(os.path.join(
            HERE, 'HR_fig2_axes_dark.png' if dark else 'HR_fig2_axes_grey.png'))
        if highlight is not None:       # fade boxes, labels and tracks toward black
            img = img.copy()
            img[..., :3] *= BG_FADE
        ax.imshow(img, extent=(X_LEFT, X_RIGHT, Y_BOT, Y_TOP),
                  aspect='auto', origin='upper', zorder=0, interpolation='bilinear')
    ax.set_xlim(X_LEFT, X_RIGHT)
    ax.set_ylim(Y_BOT, Y_TOP)

    if background:      # added stripped-He band, in the MB24 box style
        ax.add_patch(Polygon(HE_POLY, closed=True,
                             facecolor='#5c5c5c' if dark else '#cdcdcd',
                             edgecolor='none',
                             alpha=0.72*(BG_FADE if highlight is not None else 1.0),
                             zorder=1))

    def col(code):
        if highlight is not None and code not in highlight:
            return '#7a7a7a'        # everything but the highlighted codes in grey
        base = CODE_STYLE.get(code, CODE_STYLE['other'])['color']
        return DARK_COLOR.get(code, base) if dark else base

    seen = []
    if highlight is not None:       # highlighted codes drawn last, on top
        rows = sorted(rows, key=lambda r: r['code'] in highlight)
    for r in rows:
        if r['x'] is None or r.get('plot', 'yes') == 'no':
            continue
        if codes is not None and r['code'] not in codes:
            continue
        st = CODE_STYLE.get(r['code'], CODE_STYLE['other'])
        hl = highlight is None or r['code'] in highlight
        ax.plot(r['x'], r['y'], marker=st['marker'], ms=float(r.get('ms') or 11), ls='none',
                mfc=col(r['code']),
                mec=(EDGE_DARK if dark else EDGE) if hl else 'none',
                mew=1.4, alpha=1.0 if hl else 0.22, zorder=6 if hl else 5)
        if r['code'] not in seen:
            seen.append(r['code'])

    # the code names in the freed lower-left corner ARE the legend
    import matplotlib.patheffects as pe
    FIRST = {'Athena++': '2015-', 'MPI-AMRVAC': '2022-', 'AREPO': '2025-',
             'AthenaK': '2026-'}
    ytxt = 3.55
    for c in [k for k in CODE_STYLE if k in seen]:
        ax.text(5.47, ytxt, '%s (%s)' % (c, FIRST.get(c, '')), fontsize=18, color=col(c),
                alpha=1.0 if highlight is None or c in highlight else 0.6,
                ha='left', va='center', zorder=7, fontweight='bold',
                path_effects=[pe.withStroke(linewidth=2.6,
                                            foreground='black' if dark
                                            else '#2b2820')])
        ytxt -= 0.35

    ax.set_xlabel(r'$\log_{10}(T_{\rm eff}/[\mathrm{K}])$', fontsize=22, color=fg)
    ax.set_ylabel(r'$\log_{10}(L/L_\odot)$', fontsize=22, color=fg)
    ax.tick_params(axis='both', which='major', labelsize=18, colors=fg)
    if not background:
        ax.grid(alpha=0.3)

    fig.tight_layout()
    if background:      # label along the band; angle needs the final axes box
        import math
        p0 = ax.transData.transform(HE_MID[0])
        p1 = ax.transData.transform(HE_MID[1])
        ang = math.degrees(math.atan2(p1[1] - p0[1], p1[0] - p0[0]))
        if ang > 90:
            ang -= 180
        elif ang < -90:
            ang += 180
        ax.text(HE_LABEL_AT[0], HE_LABEL_AT[1], 'stripped He',
                family='serif', fontsize=25, rotation=ang,
                rotation_mode='anchor',
                color='#d0d0d0' if dark else '#8a8a8a',
                alpha=BG_FADE if highlight is not None else 1.0,
                ha='center', va='center', zorder=2)
    kw = dict(dpi=200, facecolor=fig.get_facecolor())
    if transparent:
        kw = dict(dpi=200, transparent=True)
    fig.savefig(os.path.join(HERE, out + '.png'), **kw)
    if not transparent:
        fig.savefig(os.path.join(HERE, out + '.pdf'),
                    facecolor=fig.get_facecolor())
    plt.close(fig)
    print('wrote', out + '.png/.pdf', 'points:',
          sum(1 for r in rows if r['x'] is not None), '/', len(rows))


if __name__ == '__main__':
    rows = load(os.path.join(HERE, 'sims_3d.csv'))
    # build-up sequence for slides: identical framing, growing content
    STAGES = [('hrd_3d_dark_1_background', []),
              ('hrd_3d_dark_2_codes', ['Athena++', 'MPI-AMRVAC', 'AREPO']),
              ('hrd_3d_dark_3_athenak', None)]
    for name, cs in STAGES:
        draw(rows, True, name, dark=True, codes=cs)
        draw(rows, True, name + '_transparent', dark=True, transparent=True,
             codes=cs)
    # AREPO only in colour, the other codes grey, no AthenaK points
    NOAK = ['Athena++', 'MPI-AMRVAC', 'AREPO']
    draw(rows, True, 'hrd_3d_dark_arepo', dark=True, codes=NOAK, highlight=['AREPO'])
    draw(rows, True, 'hrd_3d_dark_arepo_transparent', dark=True, transparent=True,
         codes=NOAK, highlight=['AREPO'])
    draw(rows, True, 'hrd_3d')
    draw(rows, False, 'hrd_3d_nobg')
    draw(rows, True, 'hrd_3d_dark', dark=True)
    draw(rows, True, 'hrd_3d_dark_transparent', dark=True, transparent=True)
