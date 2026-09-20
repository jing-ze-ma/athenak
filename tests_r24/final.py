"""tests_r24: corrected budget (inner face cancelled), required v_pre for the ceiling,
and the fine inner profile of the collapsed shell.  READ-ONLY."""
import re, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r24')
from state3d import load, RS, TURN, VCEIL  # noqa: E402

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
A = B + '/tests_r22/lr1'
LW = 4.0776e37
LSTAR = 2.3066e38
TAUD, TEND = 500.0, 18820.0

h = np.loadtxt(A + '/he4.hydro.hst')
th, etot, mass = h[:, 0], h[:, 6], h[:, 2]
pat = re.compile(r'erg in=(\S+) out=(\S+) \| g in=(\S+) out=(\S+)'
                 r'(?: \| L_rad,out/L = (\S+) L_rad,cut/L = (\S+))? \(t = (\S+) s\)')
rows = [[float(x) for x in (m.group(7), m.group(1), m.group(2), m.group(5) or 'nan')]
        for m in (pat.search(l) for l in open(B+'/tests_r22/lr1.log', errors='ignore'))
        if m]
tf, Ein, Eout, Lro = np.array(rows).T
ev = np.loadtxt(A + '/he4.log')
cl = np.array([(int(m.group(1)), float(m.group(2))) for m in
               (re.search(r'cycle=(\d+) time=(\S+) dt=', l)
                for l in open(B+'/tests_r22/lr1.log', errors='ignore')) if m])
tev = np.interp(ev[:, 0], cl[:, 0], cl[:, 1])
vrate = ev[:, 4]/np.maximum(np.gradient(tev), 1e-30)      # clips per second
# shell-mean damping sink from rt_profile
import struct
P = []
f = open(A+'/rt_profile.bin', 'rb')
while True:
    hh = f.read(16)
    if len(hh) < 16: break
    t, n1, nv = struct.unpack('<dii', hh)
    b = f.read(8*(n1+nv*n1))
    if len(b) < 8*(n1+nv*n1): break
    a = np.frombuffer(b, '<f8'); P.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1)))
r = P[0][1]
rf = np.empty(len(r)+1); rf[1:-1] = .5*(r[1:]+r[:-1])
rf[0] = 2*r[0]-rf[1]; rf[-1] = 2*r[-1]-rf[-2]
vol = 0.17678*4*np.pi/3*(rf[1:]**3-rf[:-1]**3)
tp = np.array([x[0] for x in P])
Pd = np.array([(P[i][2][0]*vol*P[i][2][1]**2).sum()/TAUD for i in range(len(P))])
Pd[tp >= TEND] = 0.0
# mean mass of a ceiling cell, from the dumps (g)
mbar = {3.5: 1.212e22/4928, 4.0: 7.046e21/1673}

print("CORRECTED BUDGET   (inner hydro face = 0: wall_noflux cancels it cell-by-cell;"
      " verified by the mass budget).  All in L_w = 4.0776e37 erg/s")
print("%5s %5s | %7s %7s %7s %7s %7s | %8s | %9s %7s"
      % ("t0", "t1", "dE/dt", "L_in", "-L_out", "hyd_out", "damp", "residual",
         "clips/s", "v_pre/vc"))
edges = list(np.arange(2.0, th[-1]/TURN, 0.25)) + [th[-1]/TURN]
for a, b in zip(edges[:-1], edges[1:]):
    t0, t1 = a*TURN, b*TURN
    if t1-t0 < 60: continue
    dt = t1-t0
    dE = (np.interp(t1, th, etot)-np.interp(t0, th, etot))/dt
    lo = np.nanmean(Lro[(tf >= t0) & (tf < t1)])*LSTAR
    ho = (np.interp(t1, tf, Eout)-np.interp(t0, tf, Eout))/dt
    pd = Pd[(tp >= t0) & (tp < t1)].mean() if ((tp >= t0) & (tp < t1)).any() else 0.
    res = dE - (LW - lo + ho - pd)
    cr = np.interp(0.5*(t0+t1), tev, vrate)
    mb = mbar[4.0] if b > 3.75 else mbar[3.5]
    x = -res*2.0/max(cr*mb*VCEIL**2, 1e-30)
    vp = np.sqrt(1.0+x) if x > -1 else np.nan
    print("%5.2f %5.2f | %7.3f %7.3f %7.3f %7.3f %7.3f | %8.3f | %9.2e %7.3f"
          % (a, b, dE/LW, 1.0, -lo/LW, ho/LW, -pd/LW, res/LW, cr,
             vp if res < 0 else np.nan))

print("\nFINE INNER PROFILE AT 4.0 TURNOVERS (the collapsed shell) vs 2.0")
out = {}
for idx in (4, 8):
    d, g, rc, dr, dV, rl, rr = load(idx)
    dm = (g['dens']*dV)
    out[idx] = (rc, dr, dm.sum(axis=(0, 1)), (g['dens']*dV).sum(axis=(0, 1))/dV.sum(axis=(0, 1)),
                (g['dens']*g['velx']*dV).sum(axis=(0, 1))/dm.sum(axis=(0, 1)), dm.sum())
rc = out[4][0]
print("%3s %7s | %10s %10s %9s | %10s %10s %9s %8s"
      % ("i", "r/R", "rho(2.0)", "dm(2.0)", "v1(2.0)", "rho(4.0)", "dm(4.0)", "v1(4.0)",
         "dm ratio"))
for i in list(range(0, 46, 2)):
    print("%3d %7.4f | %10.3e %10.3e %9.2e | %10.3e %10.3e %9.2e %8.2f"
          % (i, rc[i]/RS, out[4][3][i], out[4][2][i], out[4][4][i],
             out[8][3][i], out[8][2][i], out[8][4][i],
             out[8][2][i]/out[4][2][i]))
print("cumulative mass fraction inside 0.70 R: %.4f (2.0) -> %.4f (4.0)"
      % (out[4][2][rc < 0.70*RS].sum()/out[4][5], out[8][2][rc < 0.70*RS].sum()/out[8][5]))
