#!/usr/bin/env python3
"""atm_runs.py RUNDIR [TAG]: BSG atmosphere-band diagnostics on a 3-D run's bin dumps.
R_ph from the first dump with m1 data: the outermost radius where the shell-mean E >= 4 F_r/c
(Eddington tau = 2/3: E = a T_eff^4 = 4 F/c; F = L/(4 pi r^2), L the median
4 pi r^2 F_r of the inner half): the innermost radius where E <= 4 F/c; or env RPH_RSUN.  Bands (r/R_ph):
0.99-1.1, 1.01-1.1 (atmosphere), 0.5-0.9 (interior): per dump max over band radii of
|<v_r>| (mass-weighted shell mean, km/s), rms over band radii, max |<rho>/<rho>_0 - 1|.
FFT period of the shell-mean v_r at 1.0 and 1.02 R_ph (dumps after 0.5 d), when the dump
cadence dt <= 0.4 d (Nyquist period <= 0.8 d < the ~1.3 d cut-off period)."""
import sys, os, glob, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.expanduser('~/ATHENAK_vis'))
from gate_atm import shell
import bin_convert as bc
RS, DAY, C = 6.957e10, 86400.0, 2.99792458e10
run = sys.argv[1]; tag = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(run.rstrip('/'))
bd = os.path.join(run, 'bin') if os.path.isdir(os.path.join(run, 'bin')) else run
fs = sorted(glob.glob(os.path.join(bd, '*.hydro_w.*.bin')))
fm = sorted(glob.glob(os.path.join(bd, '*.m1.*.bin')))
step = int(os.environ.get('STEP', '1'))
fs = fs[::step]
rph = None
if fm:
    h = bc.read_binary_as_athdf(fm[0])
    x2f = np.asarray(h['x2f'], float)
    dom = (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :, None]*np.diff(np.asarray(h['x3f'], float))[:, None, None]
    E = (np.asarray(h['m1_e'], float)*dom).sum(axis=(0, 1))/dom.sum()
    F = (np.asarray(h['m1_f1'], float)*dom).sum(axis=(0, 1))/dom.sum()
    r = np.asarray(h['x1v'], float)
    L = np.median((4*np.pi*r**2*F)[r < r[0] + 0.5*(r[-1]-r[0])])
    ok = np.where(E <= 1.1*4*L/(4*np.pi*r**2)/C)[0]   # 1.1: an isothermal-T_ph atmosphere only touches 4F/c
    rph = r[ok[0]] if len(ok) else None
if os.environ.get('RPH_RSUN'):
    rph = float(os.environ['RPH_RSUN'])*RS
T, V, D = [], [], []
for f in fs:
    t, r, rho, vr = shell(f)
    T.append(t/DAY); V.append(vr); D.append(rho)
T = np.array(T); V = np.array(V); D = np.array(D); x = r/rph
cad = float(np.median(np.diff(T))) if len(T) > 1 else np.nan
res = dict(tag=tag, run=run, rph_rsun=rph/RS, ndump=len(T), t_end=float(T[-1]), cadence_d=cad, bands={})
print('== %s: %d dumps %.2f-%.2f d, cadence %.3f d, R_ph %.3f Rsun (%s)' % (
    tag, len(T), T[0], T[-1], cad, rph/RS, 'E=4L/(4pi r^2 c) of dump 0' if not os.environ.get('RPH_RSUN') else 'given'))
for lo, hi in [(0.99, 1.1), (1.01, 1.1), (0.5, 0.9)]:
    b = (x >= lo) & (x < hi)
    vm = np.abs(V[:, b]).max(axis=1); vr_ = np.sqrt((V[:, b]**2).mean(axis=1)); ds = np.abs(D[:, b]/D[0, b] - 1).max(axis=1)
    m1 = T >= 1.0
    res['bands']['%.2f-%.2f' % (lo, hi)] = dict(t=T.tolist(), vmax=vm.tolist(), vrms=vr_.tolist(), dsw=ds.tolist())
    print('  band %.2f-%.2f: max|<v_r>| %.2f km/s (t %.1f d), median over t>=1 d %.2f, max rms %.2f, median rms %.2f, max drho %.3g, median drho %.3g'
          % (lo, hi, vm.max(), T[int(np.argmax(vm))], np.median(vm[m1]) if m1.any() else np.nan, vr_.max(),
             np.median(vr_[m1]) if m1.any() else np.nan, ds.max(), np.median(ds[m1]) if m1.any() else np.nan))
    # mean signed v (infall indicator)
    vs = V[:, b].mean(axis=1)
    print('     signed band-mean v_r: min %.2f max %.2f km/s; fraction of dumps with mean < -1 km/s: %.2f' % (vs.min(), vs.max(), float((vs < -1).mean())))
per = {}
if cad <= 0.4:
    m = T >= 0.5
    for q in (1.0, 1.02):
        i = int(np.argmin(abs(x - q)))
        s = V[m, i] - V[m, i].mean()
        fr = np.fft.rfftfreq(m.sum(), d=cad); p = np.abs(np.fft.rfft(s*np.hanning(len(s))))**2
        k = 1 + int(np.argmax(p[1:])); per['%.2f' % q] = float(1/fr[k])
    print('  FFT periods [d]:', per)
else:
    print('  FFT: cadence %.2f d cannot resolve ~1.3 d periods (Nyquist %.2f d)' % (cad, 2*cad))
res['period_d'] = per
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '%s.json') % tag, 'w'))
