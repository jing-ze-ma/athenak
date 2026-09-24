"""ab3_ana.py <root> [ref=c2] [arm=e4]: the 3-rotation A/B of the production cadence.

<root>/<arm>/rst/*.rst every 0.5 rotation, run1.log + run2.log (ck verbose lines and the
elapsed/cycle/time diagnostics), dhj.hydro.hst (%.17e).  Per 0.5-rotation output:
  * arm vs ref: T relative difference (max, rms all/day/night, rms day/night per pressure
    band, as fast/ab_ana.py), and the same for |v| is not attempted (winds are in bin);
  * per arm: kinked columns (osc > 0.1 / > 0.3) in the 1e-5..1e-7 and < 1e-7 bar bands,
    non-converged calls, mean passes per implicit call, mean |ckdesum| in the interval;
  * energy: tot-E of the arm relative to ref at the output time, and for the cadence arm
    the interval means of edef and linerr (### ck_cadence lines).
"""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert  # noqa: E402
import dhjcs  # noqa: E402

NMB, N1, N2, N3, NG = 24, 128, 16, 16, 2
O1, O2, O3 = N1 + 2*NG, N2 + 2*NG, N3 + 2*NG
NC = O1*O2*O3
P_ROT = 3.05e5


def decode(fn):
    """hydro restart -> d, T, p on the active cells (T, p from the general-EOS caches,
    which are the first two tail slabs whatever follows them)."""
    b = open(fn, 'rb').read()
    tot = len(b)
    nv = 5
    for nx in range(2, 48):
        ds = 8*NC*(nv + nx)
        st = tot - NMB*ds
        if st < 8 or int(np.frombuffer(b[st-8:st], '<u8')[0]) != ds:
            continue
        a = np.frombuffer(b[st:], '<f8').reshape(NMB, ds//8)
        s = (slice(None), slice(NG, -NG), slice(NG, -NG), slice(NG, -NG))
        u = a[:, :nv*NC].reshape(NMB, nv, O3, O2, O1)
        tt = a[:, nv*NC:(nv+1)*NC].reshape(NMB, O3, O2, O1)
        pp = a[:, (nv+1)*NC:(nv+2)*NC].reshape(NMB, O3, O2, O1)
        return dict(d=u[:, 0][s], T=tt[s], p=pp[s])
    raise RuntimeError('layout not recognised ' + fn)


def rsttime(fn):
    d = open(fn, 'rb').read(400000)
    i = d.find(b'<par_end>')
    i = d.find(b'\n', i) + 1 + 4 + 4 + 9*8 + 2*19*4
    return float(np.frombuffer(d[i:i+8], '<f8')[0])


GEOF = '/viper/u2/jinma/ATHENAK/bench/cs_mhd_prod3/bin/dhj.mhd_w_bcc.00000.bin'
G = np.asarray(bin_convert.read_binary(GEOF)['mb_geometry'])
LON = np.empty((24, 16, 16))
for m in range(24):
    xi = np.pi/4*(G[m, 2]+(G[m, 3]-G[m, 2])*(np.arange(16)+0.5)/16)
    eta = np.pi/4*(G[m, 4]+(G[m, 5]-G[m, 4])*(np.arange(16)+0.5)/16)
    XI, ETA = np.meshgrid(xi, eta)
    X = dhjcs.panel_map(m//4, XI, ETA)
    LON[m] = np.arctan2(-X[..., 1], -X[..., 0])
DAY = np.broadcast_to((np.abs(LON) < np.pi/2)[..., None], (24, 16, 16, 128))
BANDS = [('p>1e-3', 1e30, 1e-3), ('1e-3..1e-5', 1e-3, 1e-5), ('1e-5..1e-7', 1e-5, 1e-7),
         ('<1e-7', 1e-7, 0.0)]    # bar


def kinks(T, p):
    osc = np.zeros_like(T)
    osc[..., 1:-1] = np.abs(T[..., 1:-1]-0.5*(T[..., :-2]+T[..., 2:]))/T[..., 1:-1]
    out = []
    for nm, a, b in BANDS[2:]:
        sel = (p < a*1e6) & (p >= b*1e6)
        cm = np.where(sel, osc, 0).max(-1)
        out.append((int((cm > .1).sum()), int((cm > .3).sum())))
    return out


def logs(arm):
    txt = ''
    for f in ('run1.log', 'run2.log'):
        fn = os.path.join(ROOT, arm, f)
        if os.path.exists(fn):
            txt += open(fn).read()
    cyc = [(int(c), float(t)) for c, t in re.findall(r'cycle=(\d+) time=(\S+)', txt)]
    cyc = np.array(sorted(set(cyc))) if cyc else np.zeros((0, 2))
    calls = [(int(n), int(p), abs(float(g)), 'NOT-CONV' in ln)
             for ln, n, p, g in re.findall(
                 r'(### ck_implicit ncycle=(\d+) passes=(\d+) .*?ckdesum=(\S+).*)', txt)]
    cad = [(int(n), float(e), float(le)) for n, e, le in re.findall(
        r'### ck_cadence ncycle=(\d+) .*?edef=(\S+) .*?linerr=(\S+)', txt)]
    return cyc, calls, cad


def c2t(cyc, n):
    return np.interp(n, cyc[:, 0], cyc[:, 1]) if len(cyc) else np.nan


def hst(arm):
    fs = glob.glob(os.path.join(ROOT, arm, '*.hst'))
    h = np.loadtxt(fs[0])
    return h[:, 0], h[:, 6]


ROOT = None


def main():
    global ROOT
    ROOT = sys.argv[1]
    REF = sys.argv[2] if len(sys.argv) > 2 else 'c2'
    ARM = sys.argv[3] if len(sys.argv) > 3 else 'e4'
    R = {a: {round(rsttime(f)/P_ROT, 2): f
             for f in glob.glob(os.path.join(ROOT, a, 'rst', '*.rst'))}
         for a in (REF, ARM)}
    L = {a: logs(a) for a in (REF, ARM)}
    H = {a: hst(a) for a in (REF, ARM)}
    rots = sorted(r for r in R[REF] if r in R[ARM] and abs(r*2 - round(r*2)) < 1e-6)
    print('%s - %s, T relative difference; per band rms day/night [bar]; kinks '
          '>0.1/>0.3 in 1e-5..1e-7 and <1e-7 bar per arm; nc, passes, |gap| over the '
          'interval ending there'
          % (ARM, REF))
    print('%6s %8s %8s %8s %8s' % ('rot', 'max', 'rms', 'rms day', 'rms ngt')
          + ''.join(' %17s' % b[0] for b in BANDS)
          + ' | %-15s %-15s | %9s %9s | %5s %5s %5s %5s | %8s %8s | %9s %9s %9s'
          % ('kinks ' + REF, 'kinks ' + ARM, 'dE/E', 'tE', 'nc' + REF[:2], 'nc' + ARM[:2],
             'ps' + REF[:2], 'ps' + ARM[:2], 'gap' + REF[:2], 'gap' + ARM[:2],
             'edef', 'linerr', 'lin max'))
    for r in rots:
        a = decode(R[REF][r])
        b = decode(R[ARM][r])
        d = np.abs(b['T']/a['T'] - 1.0)
        row = '%6.2f %8.1e %8.1e %8.1e %8.1e' % (r, d.max(), np.sqrt((d**2).mean()),
                                                 np.sqrt((d[DAY]**2).mean()),
                                                 np.sqrt((d[~DAY]**2).mean()))
        for nm, hi, lo in BANDS:
            sel = (a['p'] < hi*1e6) & (a['p'] >= lo*1e6)
            sd, sn = sel & DAY, sel & ~DAY
            row += ' %8.1e/%8.1e' % (np.sqrt((d[sd]**2).mean()) if sd.any() else 0,
                                     np.sqrt((d[sn]**2).mean()) if sn.any() else 0)
        ka, kb = kinks(a['T'], a['p']), kinks(b['T'], b['p'])
        row += ' | %-15s %-15s' % (' '.join('%d/%d' % k for k in ka),
                                   ' '.join('%d/%d' % k for k in kb))
        t = r*P_ROT
        ea = np.interp(t, *H[REF])
        eb = np.interp(t, *H[ARM])
        row += ' | %9.2e %9.2e' % (eb/ea - 1.0, ea)
        t0 = (r - 0.5)*P_ROT
        st = []
        for arm in (REF, ARM):
            cyc, calls, cad = L[arm]
            sel = [c for c in calls if t0 < c2t(cyc, c[0]) <= t]
            st.append((sum(c[3] for c in sel), np.mean([c[1] for c in sel]) if sel else 0,
                       np.mean([c[2] for c in sel]) if sel else 0))
        row += ' | %5d %5d %5.2f %5.2f | %8.1e %8.1e' % (st[0][0], st[1][0], st[0][1],
                                                         st[1][1], st[0][2], st[1][2])
        cyc, calls, cad = L[ARM]
        cs = [c for c in cad if t0 < c2t(cyc, c[0]) <= t]
        ed = [c[1] for c in cs if c[1] != 0.0]
        le = [c[2] for c in cs if c[2] >= 0.0]
        row += ' | %9.2e %9.2e %9.2e' % (np.mean(ed) if ed else 0,
                                         np.mean(le) if le else 0, max(le) if le else 0)
        print(row)


if __name__ == '__main__':
    main()
