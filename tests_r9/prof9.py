#!/usr/bin/env python3
"""tests_r9: the 0.97 R pile-up, off red_giant's rt_profile.bin shell means.

Usage: prof9.py <rt_profile.bin> [--turn 4705] [--shells 0.94,...] [--eint-raw]
--raw: the file was written by a binary BEFORE the tests_r9 slot-5/6 fix, so slot 5 is the
CODE temperature (kelvin/1.2027e-8) and slot 6 is eint + rho*Phi with Phi = GM(1/rin-1/r).
Both are corrected here, so an old and a new file give the same numbers.
"""
import sys, struct
import numpy as np

RSTAR = 2.3717e11
RIN   = 1.18585e11
GM    = 6.674e-8*6.2651e33
LSTAR = 2.3066e38
SIGSB = 5.670374419e-5
ACONST = 4.0*SIGSB/2.99792458e10
CLIGHT = 2.99792458e10
TURN  = 4705.0
TCGS  = 1.660538921e-24/1.3806488e-16   # EOS_Data::temp_cgs on a cgs <units> block
# the EOS radiation taper of he4_presn_cs.athinput
RHO_HI, RHO_LO = 1.93e-9, 4.97e-10
T_HI,   T_LO   = 6.3281e4, 4.5213e4
OPAC = "/viper/u2/jinma/ATHENAK/bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt"
COL  = "/viper/u2/jinma/ATHENAK/bench/hestar_presn/column_he4_presn_sph.txt"


def read_prof(fn):
    recs = []
    with open(fn, "rb") as f:
        while True:
            hd = f.read(16)
            if len(hd) < 16:
                break
            t, n1, nv = struct.unpack("<dii", hd)
            need = 8*(n1 + nv*n1)
            body = f.read(need)
            if len(body) < need:
                break
            a = np.frombuffer(body, dtype="<f8")
            recs.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1).copy()))
    return recs


def load_opac(fn=OPAC):
    hdr, vals = [], []
    for ln in open(fn):
        if ln.startswith("#"):
            hdr.append(ln)
        else:
            vals.append(float(ln))
    g = [h for h in hdr if "217" in h or h.strip().startswith("# 2")]
    nT, nD, lTmin, dlT, lDmin, dlD = [float(x) for x in g[-1].strip("# \n").split()]
    nT, nD = int(nT), int(nD)
    k = np.array(vals).reshape(nT, nD)   # T slowest
    return k, lTmin, dlT, nT, lDmin, dlD, nD


_OP = load_opac()


def kappa(T, rho):
    k, lTmin, dlT, nT, lDmin, dlD, nD = _OP
    x = (np.log10(T) - lTmin)/dlT
    y = (np.log10(rho) - lDmin)/dlD
    x = np.clip(x, 0, nT-1.000001); y = np.clip(y, 0, nD-1.000001)
    i0 = np.floor(x).astype(int); j0 = np.floor(y).astype(int)
    fx = x - i0; fy = y - j0
    lk = ((1-fx)*(1-fy)*k[i0, j0] + fx*(1-fy)*k[i0+1, j0]
          + (1-fx)*fy*k[i0, j0+1] + fx*fy*k[i0+1, j0+1])
    return 10.0**lk


def smooth(s):
    s = np.clip(s, 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def taper_w(rho, T):
    wr = smooth((np.log10(rho) - np.log10(RHO_LO))/(np.log10(RHO_HI)-np.log10(RHO_LO)))
    wt = smooth((np.log10(T) - np.log10(T_LO))/(np.log10(T_HI)-np.log10(T_LO)))
    return np.maximum(wr, wt)


def col_tau_top(rtop):
    d = np.loadtxt(COL, comments="#")
    o = np.argsort(d[:, 0])
    return float(np.interp(rtop, d[o, 0], d[o, 8]))   # col 8 = tau, r DESCENDING


def derive(r, rho, T, v1, rv1, eintc):
    n = len(r)
    kap = kappa(T, rho)
    # cell faces from the cell centres (the grid is stretched): midpoints, edges extended
    rf = np.empty(n+1)
    rf[1:-1] = 0.5*(r[:-1] + r[1:])
    rf[0] = r[0] - (rf[1]-r[0]); rf[-1] = r[-1] + (r[-1]-rf[-2])
    dr = np.diff(rf)
    # tau downward from the top, anchored on the initial column's tau at the top face
    tau = np.empty(n)
    acc = col_tau_top(rf[-1])
    for i in range(n-1, -1, -1):
        tau[i] = acc + 0.5*kap[i]*rho[i]*dr[i]
        acc += kap[i]*rho[i]*dr[i]
    # the radiative-diffusion flux on the faces, then back to centres
    lT = np.log(T); lr = np.log(r)
    dTdr = np.gradient(T, r)
    Fdiff = -(16.0*SIGSB*T**3/(3.0*kap*rho))*dTdr
    HT = np.abs(T/np.where(dTdr == 0, 1e-30, dTdr))
    F2s = Fdiff*(1.0 - HT/(2.0*r))            # the spherical-form factor (see the input)
    Freq = LSTAR/(4.0*np.pi*r**2)
    g = GM/r**2
    Gam = kap*F2s/(CLIGHT*g)
    w = taper_w(rho, T)
    Mdot = 4.0*np.pi*r**2*rv1
    return dict(r=r, rf=rf, dr=dr, rho=rho, T=T, v1=v1, Mdot=Mdot, kap=kap, tau=tau,
                Fdiff=Fdiff, F2s=F2s, Freq=Freq, Gam=Gam, w=w, eint=eintc, g=g)


def main():
    fn = sys.argv[1]
    raw = "--raw" in sys.argv          # pre-fix binary: both slots need correcting
    traw = raw or ("--traw" in sys.argv)  # slot 5 only (athena_v5, slot 6 already fixed)
    shells = [0.94, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01]
    if "--shells" in sys.argv:
        shells = [float(x) for x in sys.argv[sys.argv.index("--shells")+1].split(",")]
    recs = read_prof(fn)
    print("# %d records, t = %.1f .. %.1f s (%.3f .. %.3f turnovers)"
          % (len(recs), recs[0][0], recs[-1][0], recs[0][0]/TURN, recs[-1][0]/TURN))
    want = np.arange(0.0, 3.001, 0.25)*TURN
    picked, used = [], set()
    for tw in want:
        k = int(np.argmin([abs(rr[0]-tw) for rr in recs]))
        if k in used or abs(recs[k][0]-tw) > 0.13*TURN:
            continue
        used.add(k); picked.append(k)
    ref = None
    for k in picked:
        t, x1v, q = recs[k]
        rho, v1, rv1 = q[0], q[1], q[2]
        T = q[5]*TCGS if traw else q[5]
        eint = q[6] - (rho*GM*(1.0/RIN - 1.0/x1v) if raw else 0.0)
        D = derive(x1v, rho, T, v1, rv1, eint)
        if ref is None:
            ref = D
        print("\n=== t = %9.1f s = %.3f turnovers ===" % (t, t/TURN))
        print("  r/R      r[cm]     rho        rho/rho0   T[K]      T/T0    "
              "v1[cm/s]    Mdot[g/s]    tau      F2s/Freq  w      Gamma   eint")
        for s in shells:
            i = int(np.argmin(np.abs(D["r"]/RSTAR - s)))
            print("  %5.3f %.4e %.4e %8.3f  %.4e %7.3f %+.3e %+.3e %.3e %8.3f %6.4f %7.3f %.3e"
                  % (D["r"][i]/RSTAR, D["r"][i], D["rho"][i], D["rho"][i]/ref["rho"][i],
                     D["T"][i], D["T"][i]/ref["T"][i], D["v1"][i], D["Mdot"][i],
                     D["tau"][i], D["F2s"][i]/D["Freq"][i], D["w"][i], D["Gam"][i],
                     D["eint"][i]))
    # the shell budget of the 0.97 R cell, between its own two faces
    print("\n=== SHELL BUDGET of the 0.97 R cell (between its own faces) ===")
    print("# t/turn  rf_in       rf_out     dM/dt_in-out[g/s]  radiative div"
          "[erg/cm3/s]  PdV[erg/cm3/s]  d eint/dt_implied   eint")
    prev = None
    for k in picked:
        t, x1v, q = recs[k]
        rho, v1, rv1 = q[0], q[1], q[2]
        T = q[5]*TCGS if traw else q[5]
        eint = q[6] - (rho*GM*(1.0/RIN - 1.0/x1v) if raw else 0.0)
        D = derive(x1v, rho, T, v1, rv1, eint)
        i = int(np.argmin(np.abs(D["r"]/RSTAR - 0.97)))
        rf, dr = D["rf"], D["dr"]
        A = 4.0*np.pi*rf**2
        V = (4.0*np.pi/3.0)*(rf[i+1]**3 - rf[i]**3)
        f2sf = np.interp(rf, D["r"], D["F2s"])
        raddiv = -(A[i+1]*f2sf[i+1] - A[i]*f2sf[i])/V
        vf = np.interp(rf, D["r"], D["v1"])
        # p = gas+tapered-rad pressure is not in the record; use the divergence of v only
        divv = (A[i+1]*vf[i+1] - A[i]*vf[i])/V
        mdf = np.interp(rf, D["r"], D["Mdot"])
        print("  %6.3f  %.4e %.4e  %+.4e        %+.4e      divv=%+.3e   %.4e"
              % (t/TURN, rf[i], rf[i+1], mdf[i]-mdf[i+1], raddiv, divv, D["eint"][i]))
        if prev is not None:
            dt = t - prev[0]
            print("        measured d eint/dt = %+.4e erg/cm3/s over dt = %.1f s,"
                  "  d rho/dt = %+.4e" % ((D["eint"][i]-prev[1])/dt, dt,
                                          (D["rho"][i]-prev[2])/dt))
        prev = (t, D["eint"][i], D["rho"][i])


if __name__ == "__main__":
    main()
