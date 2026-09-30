"""IC profile from a ck RCE (rce_w121.py output), built as x1/ic_w121_1x.txt was:
RCE rows in ascending p, isothermal rows (same spacing) up to <= 1e-10 bar, and below
the RCE bottom the exact EOS adiabat (RK2 in ln p, nabla_ad from the EOS dump) to >= 2000 bar.
usage: make_ic.py <rce.txt> <eos dump dir> <out> <label>
       make_ic.py --regress   (rebuilds the 1x IC and compares with x1/ic_w121_1x.txt)"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dclib  # noqa: E402


def build(rce, dumpdir):
    r = np.loadtxt(rce)
    hdr = open(rce).readline().strip()
    lp, T = r[::-1, 0], r[::-1, 1]
    d = np.mean(np.diff(lp))
    eos = dclib.EOS(dumpdir + '/eos_grid.txt')
    top = []
    x = lp[0]
    while x > -4.0 + 1e-9:
        x -= d
        top.append(x)
    top = np.array(top[::-1])
    bot, Tb = [], []
    x, t = lp[-1], T[-1]
    h = d*np.log(10.0)

    def gad(lpb, tt):
        p = np.array([10**lpb])
        tt = np.array([tt])
        return float(eos.at('gad', eos.rho_pt(p, tt), tt)[0])
    while x < np.log10(2000.0e6) - 1e-9:
        k1 = gad(x, t)
        tm = t*np.exp(0.5*h*k1)
        k2 = gad(x + 0.5*d, tm)
        t = t*np.exp(h*k2)
        x += d
        bot.append(x)
        Tb.append(t)
    lpo = np.concatenate([top, lp, bot])
    To = np.concatenate([np.full(len(top), T[0]), T, Tb])
    return lpo, To, hdr, d, lp[-1]


if sys.argv[1] == '--regress':
    X = '/viper/ptmp2/jinma/wasp121_0925/x1/'
    lpo, To, hdr, d, lb = build(X + 'rce/ckrce_nq2_B.txt', '/viper/ptmp2/jinma/deepconv_0925/dump')
    ref = np.loadtxt(X + 'ic_w121_1x.txt')
    print('rows', len(lpo), len(ref), 'max |d lp|', np.abs(lpo - ref[:, 0]).max(),
          'max |dT/T|', np.abs(To/ref[:, 1] - 1).max())
    sys.exit()
rce, dumpdir, out, label = sys.argv[1:5]
lpo, To, hdr, d, lb = build(rce, dumpdir)
with open(out, 'w') as fh:
    fh.write('# WASP-121b IC, %s: 1-D global-mean ck RCE, route B, ck_nquad 2, 20 cells/dex, '
             'Bond albedo 0.277 (w121prod_0929/met3, make_ic.py)\n' % label)
    fh.write('# UNIFORM rows (%.8f dex): isothermal to 1e-10 bar at the top; below the RCE bottom '
             '(%.0f bar) the exact EOS adiabat (RK2) to >= 2000 bar.\n' % (d, 10**lb/1e6))
    fh.write('# solver header: %s\n' % hdr.lstrip('# '))
    fh.write('# log10 p [barye]  T [K]\n')
    for a, b in zip(lpo, To):
        fh.write('%.8f %.6f\n' % (a, b))
print('rows', len(lpo), 'bottom', 10**lpo[-1]/1e6, 'bar', To[-1], 'K')
