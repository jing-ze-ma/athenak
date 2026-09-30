#!/usr/bin/env python3
"""oe.py: radial odd-even (2-dr) diagnostics per restart, same definitions as
../deepmix_c32/deepmix.py (oe = shell rms(d2 v_r)/4 / shell rms(v_r), dOmega-weighted;
1 = pure 2-cell radial mode, ~0.61 = white noise).
usage: oe.py LABEL rst1 [rst2 ...]   -> prints one block per restart, appends to oe_<LABEL>.txt
Per restart, per pressure band: mean of per-shell oe, rms v_r [m/s], Mdot sign-alternation
fraction (share of adjacent cell shells whose l=0 Mdot changes sign), l=0 oe of Mdot
(rms d2 Mdot/4 / rms Mdot), and the median lhllc phi = chi(2-chi), chi = 2|v_r|/c,
c = sqrt(1.4 p/rho) (approximate; the solver uses max(|vl|,|vr|)/max(cl,cr)).
"""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhjcs        # noqa: E402
import dhj_remap as R  # noqa: E402

BIN = '/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.00150.bin'  # horizontal geometry
BANDS = [(100, 1e4), (30, 100), (10, 30), (1, 10), (0.1, 1), (0.01, 0.1)]
PR = 1.101535e5


def dom_of(nmb, nk, nj):
    raw = dhjcs.bin_convert.read_binary(BIN)
    Gm = np.asarray(raw['mb_geometry'])
    dom = np.empty((nmb, nk, nj))
    for m in range(nmb):
        x = np.tan(np.pi/4*(Gm[m, 2] + (Gm[m, 3] - Gm[m, 2])*(np.arange(nj) + 0.5)/nj))
        y = np.tan(np.pi/4*(Gm[m, 4] + (Gm[m, 5] - Gm[m, 4])*(np.arange(nk) + 0.5)/nk))
        X, Y = np.meshgrid(x, y)
        dxi = np.pi/4*(Gm[m, 3] - Gm[m, 2])/nj
        deta = np.pi/4*(Gm[m, 5] - Gm[m, 4])/nk
        dom[m] = dxi*deta*(1 + X*X)*(1 + Y*Y)/(1 + X*X + Y*Y)**1.5
    return dom


def main():
    lab = sys.argv[1]
    dom = None
    out = open('oe_%s.txt' % lab, 'a')
    for fn in sys.argv[2:]:
        h = R.read_rst(fn)
        ng = h['ng']
        s = (slice(None), slice(ng, -ng), slice(ng, -ng), slice(ng, -ng))
        rho = h['u'][:, 0][s]
        m1 = h['u'][:, 1][s]
        p = h['p'][s]
        if dom is None:
            dom = dom_of(*rho.shape[:3])
            wn = (dom/dom.sum())[..., None]
            nx1, r0, r1, c = R.grid_of(h['par'])
            rc = R.centroids(R.edges(nx1, r0, r1, c))
        vr = m1/rho
        vr2 = (vr*vr*wn).sum((0, 1, 2))
        d2 = np.full(vr2.shape, np.nan)
        d2[1:-1] = ((vr[..., 2:] - 2*vr[..., 1:-1] + vr[..., :-2])**2*wn).sum((0, 1, 2))/16
        oe = np.sqrt(d2/vr2)
        pb = (p*wn).sum((0, 1, 2))/1e6
        cs = np.sqrt(1.4*p/rho)
        chi = np.minimum(1.0, 2*np.abs(vr)/cs)
        phi = np.median((chi*(2 - chi)).reshape(-1, vr.shape[-1]), axis=0)
        md = rc**2*(dom[..., None]*m1).sum((0, 1, 2))
        rot = h['t']/PR
        line = '%s %s rot %.3f |' % (lab, fn.split('/')[-1], rot)
        for lo, hi in BANDS:
            k = np.where((pb >= lo) & (pb < hi) & np.isfinite(oe))[0]
            if len(k) < 3:
                line += ' %g-%g: n/a |' % (lo, hi)
                continue
            mk = md[k[0]:k[-1] + 1]
            alt = np.mean(np.sign(mk[1:]) != np.sign(mk[:-1]))
            d2m = mk[2:] - 2*mk[1:-1] + mk[:-2]
            oem = np.sqrt(np.mean(d2m**2)/16/np.mean(mk**2))
            line += ' %g-%g bar: oe %.2f vr %.2f alt %.2f oeM %.2f phi %.1e |' % (
                lo, hi, oe[k].mean(), np.sqrt(vr2[k].mean())/100, alt, oem, np.median(phi[k]))
        print(line, flush=True)
        out.write(line + '\n')


if __name__ == '__main__':
    main()
