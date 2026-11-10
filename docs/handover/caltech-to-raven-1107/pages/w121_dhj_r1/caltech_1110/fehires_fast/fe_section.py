"""fe_section.py: HTML snippet for the dhj page from the fe_series.py outputs.

  python fe_section.py [--series DIR] [--img-prefix img/] > fe_section.html
Uses <series>/fe_series.json and the images that `fe_series.py plot` wrote (fe_vs_rot.png and
fehires_vs_seidel.png; the page refresh copies them into the page img/ dir).
"""
import argparse
import html
import json
import os
import sys
import numpy as np

FF = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, FF)
import fe_series  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--series', default=fe_series.SER)
    ap.add_argument('--img-prefix', default='img/')
    a = ap.parse_args()
    S = json.load(open(os.path.join(a.series, 'fe_series.json')))
    win, ver = [], []
    for k, s in S.items():
        if not s['rot']:
            continue
        o = np.argsort(s['rot'])
        b6 = np.array(s['pc6'])[o][:, 1]
        r = np.array(s['rot'])[o]
        w = b6[-5:].mean()
        win.append('%.1f (%s%s)' % (w, s['name'], '' if s['final'] else ', rot %.0f' % r[-1]))
        v, sl = fe_series.verdict(s)
        ver.append('%s %s%s' % (s['name'], v, '' if np.isnan(sl) else ' (%+.2f km/s per 10 rot)' % sl))
    obs = ' / '.join('%.2f' % s['seidel'] for s in S.values() if s['rot'])
    p = a.img_prefix
    out = [
        '<h2 id="fehires">Fe I terminator wind vs Seidel et al. 2026</h2>',
        '<p><strong>Model Fe I net blueshift (Kp-marginalised, Fe cut above 1e-6 bar, mean of the newest 5 dumps): '
        + html.escape(', '.join(win)) + ' km/s vs observed ' + obs + ' km/s.</strong> Series: '
        + html.escape('; '.join(ver)) + '.</p>',
        '<figure><img src="%sfehires_vs_seidel.png" alt="Fe I net blueshift vs T_eq: model vs Seidel 2026" '
        'loading="lazy"><figcaption>Net Fe I blueshift (km/s, positive = blueshift) vs T_eq. Points with error '
        'bars: Seidel et al. 2026 (ESPRESSO/MAROON-X). Model: fe_hires (high-resolution Fe I transmission through '
        'the 3-D dumps, cross-correlation, Kp-marginalised offset), mean of the newest 5 dumps; filled = Fe removed '
        'above 1e-6 bar, open = no top cut; bars = 16-84 %% range over Kp.</figcaption></figure>' % p,
        '<figure><img src="%sfe_vs_rot.png" alt="Fe I blueshift vs rotation per planet" loading="lazy">'
        '<figcaption>Kp-marginalised Fe I blueshift vs rotation (one dump every 10 rotations plus the newest 5). '
        'Filled = Fe cut above 1e-6 bar, open = no cut; grey band = Seidel+2026 1 sigma. "converged": the '
        'least-squares slope over the last 30 rotations changes the value by less than max(0.2 km/s, 2 x the '
        'dump-to-dump sd) over 30 rotations; otherwise "drifting".</figcaption></figure>' % p]
    print('\n'.join(out))


if __name__ == '__main__':
    main()
