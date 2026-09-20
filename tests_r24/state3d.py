"""tests_r24: 3-D state + velocity-ceiling census for tests_r22/lr1.  READ-ONLY."""
import os, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc  # noqa: E402

RUN = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r22/lr1'
RS = 2.3717e11
TURN = 4705.0
CPOLY = [0.992525, -0.404821, -0.372255, -1.499759]
ARAD = 7.5657332503e-15
KB = 1.380649e-16
MH = 1.6605390e-24
MU = 4.0026/3.0
VCEIL = 1.0e8
RAD_LO, RAD_HI = -11.0, -10.0     # OV: eos_rad_rho_lo=1e-11, hi=1e-10
KES = 0.2                          # cm^2/g, He electron scattering (X=0)


def stretch(xi):
    u = xi.copy()
    xik = xi.copy()
    for kk in range(1, 5):
        u = u + CPOLY[kk-1]*xik*(1.0-xi)
        xik = xik*xi
    return u


def load(idx):
    d = bc.read_binary(os.path.join(RUN, 'bin', 'he4.hydro_w.%05d.bin' % idx))
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    out = {}
    nj, nk = d['nx2_mb'], d['nx3_mb']
    for nm in d['var_names']:
        g = np.empty((n3, n2, n1), dtype=np.float64)
        src = d['mb_data'][nm]
        for m in range(d['n_mbs']):
            _, lx2, lx3, _ = d['mb_logical'][m]
            g[lx3*nk:(lx3+1)*nk, lx2*nj:(lx2+1)*nj, :] = src[m]
        out[nm] = g
    rf = d['x1min'] + (d['x1max']-d['x1min'])*stretch(np.linspace(0.0, 1.0, n1+1))
    rl, rr = rf[:-1], rf[1:]
    q = rl/rr
    rc = 0.25*(q*q+1.0)/((1./3.)*(q*q+q+1.))*(rr+rl)
    dr = rr-rl
    thf = np.linspace(d['x2min'], d['x2max'], n2+1)
    wj = np.cos(thf[:-1]) - np.cos(thf[1:])          # solid angle / dphi per theta row
    dphi = (d['x3max']-d['x3min'])/n3
    dV = ((rr**3-rl**3)/3.0)[None, None, :]*wj[None, :, None]*dphi   # (1,n2,n1)
    dV = np.broadcast_to(dV, (n3, n2, n1))
    return d, out, rc, dr, dV, rl, rr


def weight(lrho):
    s = (lrho - RAD_LO)/(RAD_HI - RAD_LO)
    return np.where(s <= 0, 0.0, np.where(s >= 1, 1.0, s*s*(3-2*s)))


def report(idx):
    d, g, rc, dr, dV, rl, rr = load(idx)
    t = d['time']
    rho, vx, vy, vz, ei = g['dens'], g['velx'], g['vely'], g['velz'], g['eint']
    dm = rho*dV
    Mtot = dm.sum()
    v2 = vx*vx + vy*vy + vz*vz
    vmag = np.sqrt(v2)
    print("\n================ dump %d   t = %.1f s = %.3f turnovers" % (idx, t, t/TURN))
    print("  total mass in wedge = %.6e g   KE = %.4e erg  eint = %.4e erg"
          % (Mtot, (0.5*dm*v2).sum(), (ei*dV).sum()))
    # ---- velocity-ceiling census
    for thr in (0.5, 0.9, 0.99, 0.999999):
        sel = vmag > thr*VCEIL
        n = int(sel.sum())
        kem = (0.5*dm*v2)[sel].sum()
        print("  |v| > %8.4f vceil : %8d cells (%.3e of all)  KE = %.4e erg"
              " mass = %.3e g" % (thr, n, n/vmag.size, kem, dm[sel].sum()))
    sel = vmag > 0.99*VCEIL
    if sel.any():
        ii = np.where(sel)[2]
        print("     their r/R: min %.3f  p10 %.3f  median %.3f  p90 %.3f  max %.3f"
              % tuple(rc[np.percentile(ii, p).astype(int)]/RS
                      for p in (0, 10, 50, 90, 100)))
        print("     their rho: median %.3e  max %.3e ; KE at ceiling"
              " 0.5 m vceil^2 = %.4e erg"
              % (np.median(rho[sel]), rho[sel].max(), 0.5*dm[sel].sum()*VCEIL**2))
        # radial histogram of the count
        h = np.bincount(ii, minlength=len(rc))
        top = np.argsort(h)[-6:][::-1]
        print("     top radial bins (i, r/R, count, <rho>):",
              ", ".join("(%d,%.3f,%d,%.1e)" % (i, rc[i]/RS, h[i],
                        np.median(rho[:, :, i][vmag[:, :, i] > 0.99*VCEIL]))
                        for i in top if h[i] > 0))
    # ---- shell means
    W = dV.sum(axis=(0, 1))
    rho_s = (rho*dV).sum(axis=(0, 1))/W
    v1_s = (rho*vx*dV).sum(axis=(0, 1))/(rho*dV).sum(axis=(0, 1))
    ei_s = (ei*dV).sum(axis=(0, 1))/W
    menc = np.cumsum(dm.sum(axis=(0, 1)))
    # temperature from the shell-mean e (gas + tapered radiation), ideal-ionized He
    lr = np.log10(np.maximum(rho_s, 1e-300))
    w = weight(lr)
    lo, hi = np.full(len(rc), 3.0), np.full(len(rc), 7.0)
    for _ in range(60):
        mid = 0.5*(lo+hi)
        T = 10**mid
        e = 1.5*rho_s*KB*T/(MU*MH) + w*ARAD*T**4
        lo = np.where(e < ei_s, mid, lo)
        hi = np.where(e < ei_s, hi, mid)
    T = 10**(0.5*(lo+hi))
    pg = rho_s*KB*T/(MU*MH)
    pr = w*ARAD*T**4/3.0
    beta = pg/(pg+pr)
    tau = np.cumsum((rho_s*dr)[::-1])[::-1]*KES
    return dict(t=t, r=rc, dr=dr, rho=rho_s, v1=v1_s, T=T, beta=beta, menc=menc,
                Mtot=Mtot, tau=tau, mass_i=dm.sum(axis=(0, 1)))


if __name__ == '__main__':
    idxs = [int(a) for a in sys.argv[1:]] or [4, 6, 8]
    out = [report(i) for i in idxs]
    print("\n---- radial profiles (shell means; tau from kappa_es = 0.2) ----")
    hdr = "  i  r/R    dr/R  "
    for o in out:
        hdr += " | %5.2ft: rho      T        v1       beta   M(<r)/M  tau" % (o['t']/TURN)
    print(hdr)
    r = out[0]['r']
    sel = list(range(0, len(r), 6)) + [len(r)-1]
    for i in sorted(set(sel)):
        s = "%3d %6.4f %6.4f" % (i, r[i]/RS, out[0]['dr'][i]/RS)
        for o in out:
            s += " | %8.2e %8.2e %8.1e %6.3f %7.4f %8.2e" % (
                o['rho'][i], o['T'][i], o['v1'][i], o['beta'][i],
                o['menc'][i]/o['Mtot'], o['tau'][i])
        print(s)
    print("\n---- where the photosphere is (tau_es = 1) and the mass shifts ----")
    for o in out:
        k = int(np.argmin(np.abs(np.log(o['tau']/1.0))))
        print("  t=%.2f turn: tau=1 at i=%d r/R=%.4f rho=%.2e T=%.3e | "
              "M(wedge)=%.6e | M(r<0.55R)=%.4e M(r>1.0R)=%.4e"
              % (o['t']/TURN, k, o['r'][k]/RS, o['rho'][k], o['T'][k], o['Mtot'],
                 o['mass_i'][o['r'] < 0.55*RS].sum(),
                 o['mass_i'][o['r'] > 1.0*RS].sum()))
