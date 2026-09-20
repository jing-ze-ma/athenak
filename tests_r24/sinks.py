"""tests_r24: mass budget, energy decomposition, sponge + vceil sink sizing.  READ-ONLY."""
import os, re, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r24')
from state3d import load, weight, RS, TURN, VCEIL, ARAD, KB, MH, MU  # noqa: E402

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
A = B + '/tests_r22/lr1'
LW = 4.0776e37
LSTAR = 2.3066e38
GM = 6.674e-8*6.2651e33
RIN = 1.18585e11
SPZ = 0.96
SPC = 0.1

# ---------------- mass budget ----------------
h = np.loadtxt(A + '/he4.hydro.hst')
th, mass = h[:, 0], h[:, 2]
pat = re.compile(r'erg in=(\S+) out=(\S+) \| g in=(\S+) out=(\S+)'
                 r'(?: \| L_rad,out/L = (\S+) L_rad,cut/L = (\S+))? \(t = (\S+) s\)')
rows = [[float(x) for x in (m.group(7), m.group(1), m.group(2), m.group(3), m.group(4),
                            m.group(5) or 'nan', m.group(6) or 'nan')]
        for m in (pat.search(l) for l in open(B+'/tests_r22/lr1.log', errors='ignore'))
        if m]
F = np.array(rows)
tf, Ein, Eout, Min, Mout, Lro, Lrc = F.T
print("MASS BUDGET (g/s)   [inner face flux is cancelled cell-by-cell by wall_noflux]")
print("%6s %6s | %10s %10s %10s %10s" % ("t0", "t1", "dM/dt", "face_in", "face_out",
                                         "unaccounted"))
for a in np.arange(2.0, 4.30, 0.25):
    b = min(a+0.25, th[-1]/TURN)
    if b <= a: break
    t0, t1 = a*TURN, b*TURN
    dt = t1-t0
    dM = (np.interp(t1, th, mass)-np.interp(t0, th, mass))/dt
    mi = (np.interp(t1, tf, Min)-np.interp(t0, tf, Min))/dt
    mo = (np.interp(t1, tf, Mout)-np.interp(t0, tf, Mout))/dt
    print("%6.2f %6.2f | %10.2e %10.2e %10.2e %10.2e" % (a, b, dM, mi, mo, dM-mo))

# ---------------- event counters vs time ----------------
ev = np.loadtxt(A + '/he4.log')
cyc = ev[:, 0]
# cycle -> time from the ndiag lines
cl = [(int(m.group(1)), float(m.group(2))) for m in
      (re.search(r'cycle=(\d+) time=(\S+) dt=', l)
       for l in open(B+'/tests_r22/lr1.log', errors='ignore')) if m]
cl = np.array(cl)
tev = np.interp(cyc, cl[:, 0], cl[:, 1])
dcyc = np.gradient(cyc)
dtev = np.gradient(tev)
names = ['cycle', 'dfloor', 'efloor', 'tfloor', 'vceil', 'fail', 'c2p_it', 'fofc',
         'efloor_de', 'tclamp', 'tset']
print("\nEVENT COUNTERS (per interval) and rates")
print("%8s %7s | %10s %10s %10s %10s %10s" %
      ("turn", "dt[s]", "vceil/s", "vceil/stage", "dfloor/s", "fofc/s",
       "efloor_de[erg/s]"))
for k in range(0, len(cyc), max(1, len(cyc)//24)):
    if dtev[k] <= 0: continue
    print("%8.3f %7.2f | %10.3e %10.3e %10.3e %10.3e %10.3e" %
          (tev[k]/TURN, dtev[k]/dcyc[k], ev[k, 4]/dtev[k], ev[k, 4]/(2*dcyc[k]),
           ev[k, 1]/dtev[k], ev[k, 7]/dtev[k], ev[k, 8]/dtev[k]))

# ---------------- per-dump energy decomposition + sponge power ----------------
print("\nENERGY DECOMPOSITION FROM THE DUMPS (wedge integrals, erg)")
print("%6s %12s %12s %12s %12s %12s" % ("turn", "eint", "KE", "PE=+rho Phi",
                                        "sum=tot-E", "hst tot-E"))
prev = None
for idx in (4, 5, 6, 7, 8):
    d, g, rc, dr, dV, rl, rr = load(idx)
    t = d['time']
    rho, vx, vy, vz, ei = g['dens'], g['velx'], g['vely'], g['velz'], g['eint']
    dm = rho*dV
    v2 = vx*vx+vy*vy+vz*vz
    Ei = (ei*dV).sum()
    Ek = (0.5*dm*v2).sum()
    phi = GM*(1.0/RIN - 1.0/rc)
    Ep = (dm*phi[None, None, :]).sum()
    print("%6.3f %12.5e %12.5e %12.5e %12.5e %12.5e" %
          (t/TURN, Ei, Ek, Ep, Ei+Ek+Ep, np.interp(t, th, h[:, 6])))
    if prev is not None:
        dt = t-prev[0]
        print("        d/dt in L_w:  eint %+7.3f   KE %+7.3f   PE %+7.3f   tot %+7.3f"
              % ((Ei-prev[1])/dt/LW, (Ek-prev[2])/dt/LW, (Ep-prev[3])/dt/LW,
                 (Ei+Ek+Ep-prev[1]-prev[2]-prev[3])/dt/LW))
    prev = (t, Ei, Ek, Ep)
    # --- sponge power (radius sponge, r > x1min + 0.96*(x1max-x1min))
    zs = d['x1min'] + SPZ*(d['x1max']-d['x1min'])
    lr = np.log10(np.maximum(rho, 1e-300))
    w = weight(lr)
    # T from eint (same inversion as state3d, vectorised, coarse but enough for cs)
    lo, hi = np.full(rho.shape, 3.0), np.full(rho.shape, 7.0)
    for _ in range(40):
        mid = 0.5*(lo+hi)
        T = 10**mid
        e = 1.5*rho*KB*T/(MU*MH) + w*ARAD*T**4
        lo = np.where(e < ei, mid, lo)
        hi = np.where(e < ei, hi, mid)
    T = 10**(0.5*(lo+hi))
    pg = rho*KB*T/(MU*MH)
    pr = w*ARAD*T**4/3.0
    p = pg+pr
    beta = pg/p
    gam1 = beta*(5./3.) + (1-beta)*(4./3.)
    cs = np.sqrt(gam1*p/np.maximum(rho, 1e-300))
    ramp = np.where(rc > zs, ((rc-zs)/(d['x1max']-zs))**2, 0.0)[None, None, :]
    rate = SPC*ramp*cs/dr[None, None, :]           # 1/s;  fac = 1/(1+rate*bdt)
    Pspg = (2.0*rate*0.5*dm*v2).sum()              # (1-fac^2)/bdt -> 2*rate for small
    msp = dm[np.broadcast_to(rc > zs, rho.shape)].sum()
    print("        sponge: zs/R=%.4f  mass above %.3e g  KE above %.3e erg"
          "  P_sponge = %.3e erg/s = %.4f L_w"
          % (zs/RS, msp, (0.5*dm*v2)[np.broadcast_to(rc > zs, rho.shape)].sum(),
             Pspg, Pspg/LW))
    # --- vceil census-based sink: N at ceiling and the band population
    vm = np.sqrt(v2)
    nceil = int((vm > 0.999999*VCEIL).sum())
    Mceil = dm[vm > 0.999999*VCEIL].sum()
    band = (vm > 0.99*VCEIL) & (vm <= 0.999999*VCEIL)
    nband = int(band.sum())
    if nceil:
        dv = nceil/max(nband, 1)*0.01*VCEIL     # drift per stage from the flux balance
        P = Mceil*VCEIL*dv*2.0/np.interp(t, th, h[:, 1])
        print("        vceil: N_at_ceiling=%d  M=%.3e g  N in (0.99,1)=%d"
              "  -> dv/stage=%.2e cm/s  P_vceil~%.3e = %.3f L_w  (KE_ceil=%.3e)"
              % (nceil, Mceil, nband, dv, P, P/LW, 0.5*Mceil*VCEIL**2))
        print("        vceil UPPER bound (whole ceiling KE removed each stage):"
              " %.3f L_w" % (2*0.5*Mceil*VCEIL**2/np.interp(t, th, h[:, 1])/LW))
