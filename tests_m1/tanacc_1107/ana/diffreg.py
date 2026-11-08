"""On the cells where the tensor key changes the answer (|E_on - E_off| > thr J0): errors of each arm against
the exact solution.  python3 diffreg.py <on_run> [thr]"""
import sys, os, math, numpy as np
from evalrun import load, cfg
import exact
on = sys.argv[1].rstrip('/'); off = on.replace('_on', '_off'); thr = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
p = cfg(on); a, b, t = load(on, 'm1'), load(off, 'm1'), load(on, 'm1_lat'); tb = load(off, 'm1_lat')
r = t['r']; rh = np.stack([t['m1_rx'], t['m1_ry'], t['m1_rz']], 1)
J0 = 2*(1 - np.sqrt(1 - 1/r**2))
idx = np.where(np.abs(a['m1_e'] - b['m1_e']) > thr*J0)[0]
th, ph, rc = float(p['clump_th']), float(p['clump_ph']), float(p['clump_r'])
xc = rc*np.array([math.sin(th)*math.cos(ph), math.sin(th)*math.sin(ph), math.cos(th)])
B = float(p['arad'])*float(p['clump_temp'])**4
cl = (xc, float(p['clump_a']), 1e-30*float(p['clump_rho']), B)
X = r[idx, None]*rh[idx]
J, H, K = exact.moments(X, cl)
Hr = np.einsum('pi,pi->p', H, rh[idx])
dc = np.linalg.norm(X - xc, axis=1)/cl[1]
tau_in = np.where(dc < 1, 1, 0)
for nm, sel in (('all', dc > -1), ('inside clump', dc < 1), ('shell 1-1.5 a', (dc >= 1) & (dc < 1.5)),
                ('outside 1.5 a', dc >= 1.5)):
    if sel.sum() == 0: continue
    w = np.maximum(J0[idx][sel], J[sel])
    eon = np.abs(a['m1_e'][idx][sel] - J[sel])/w; eof = np.abs(b['m1_e'][idx][sel] - J[sel])/w
    fon = np.abs(a['m1_f1'][idx][sel] - Hr[sel])*r[idx][sel]**2; fof = np.abs(b['m1_f1'][idx][sel] - Hr[sel])*r[idx][sel]**2
    print('%-26s %-14s N %5d | E err on %.3e off %.3e (on closer in %.0f %%) | r2 Fr err on %.3e off %.3e' % (
        os.path.basename(on), nm, sel.sum(), eon.mean(), eof.mean(), 100*np.mean(eon < eof), fon.mean(), fof.mean()))
