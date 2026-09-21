"""Gates (d) and (e) for problem/ck_beam_sph, from the correlated-k column dumps.

Dump columns: 0 i  1 r_face  2 p[bar]  3 T[K]  4 F_lw_net  5 Q_sw  6 Gamma_1  7 grad_ad
              8 tau_R  9 w_diff  10 Src_lw  11 A_cell  12 V_cell
"""
import sys
import numpy as np

SIG = 5.6704e-5
TIRR = 2500.0*np.sqrt(2.0)
FSTAR = SIG*TIRR**4


def load(path):
    hdr = {}
    with open(path) as f:
        for ln in f:
            if not ln.startswith("#"):
                break
            for k in ["mu0", "icut"]:
                if k + " =" in ln:
                    hdr[k] = float(ln.split(k + " =")[1].split(",")[0]
                                   .split("(")[0].strip())
    return hdr, np.loadtxt(path)


def col(tag):
    out = {}
    for w in ["old", "new"]:
        hdr, d = load(f"bd_{tag}_{w}/col.txt")
        i = d[:, 0].astype(int)
        rf, p, Q, Ac, V = d[:, 1], d[:, 2], d[:, 5], d[:, 11], d[:, 12]
        lo = int(np.where(i == int(hdr["icut"]))[0][0])
        hi = len(i) - 1
        dz = V/np.where(Ac > 0, Ac, 1.0)
        # absorbed power per unit area of the TOP face, i.e. the column integral of the
        # volumetric rate weighted by the shell thickness -- the quantity that is
        # F* mu0 (1 - e^-tau) in the plane-parallel limit
        colabs = float(np.sum(Q[lo:hi]*dz[lo:hi]))
        nz = np.where(Q[lo:hi] > 1.0e-6*max(Q[lo:hi].max(), 1.0e-30))[0]
        out[w] = dict(mu0=hdr["mu0"], p=p, Q=Q, dz=dz, lo=lo, hi=hi, abs=colabs,
                      ptop=(p[lo + nz[0]] if nz.size else np.nan),
                      pbot=(p[lo + nz[-1]] if nz.size else np.nan),
                      ppk=(p[lo + int(np.argmax(Q[lo:hi]))] if nz.size else np.nan))
    return out


TAGS = [("p979", "+0.979"), ("p435", "+0.435"), ("p102", "+0.102"),
        ("n102", "-0.102"), ("n195", "-0.195"), ("n513", "-0.513")]

if __name__ == "__main__":
    print("(d) HEATING PROFILES, production radial grid, t = 0")
    print("     mu0    | absorbed per unit top area [erg/s/cm2]  | "
          "/(F* |mu0|)     | peak-heating p [bar]  | heating spans p [bar]")
    rows = {}
    for tag, lab in TAGS:
        c = col(tag)
        rows[tag] = c
        o, n = c["old"], c["new"]
        mu = abs(o["mu0"])
        print(f"  {lab:>7} | old {o['abs']:11.4e}  new {n['abs']:11.4e} | "
              f"old {o['abs']/(FSTAR*mu):7.4f} new {n['abs']/(FSTAR*mu):7.4f} | "
              f"old {o['ppk']:9.3e} new {n['ppk']:9.3e} | "
              f"new {n['ptop']:8.2e}..{n['pbot']:8.2e}")
    print(f"  F* = {FSTAR:.5e} erg/s/cm2,  (1-albedo) = 0.991541")
    print("\n(e) GLOBAL ABSORBED POWER over the sphere.")
    print("    P = int over the lit hemisphere of (absorbed per unit top area) dA.")
    print("    The old scheme absorbs only where mu0 > 0 and clamps 1/mu0 at 10;")
    print("    the new one also lights the twilight ring.  Using the measured")
    print("    column integrals as a function of mu0 and integrating 2 pi r^2 "
          "d(mu0):")
    mus, ao, an = [], [], []
    for tag, lab in TAGS:
        c = rows[tag]
        mus.append(c["old"]["mu0"])
        ao.append(c["old"]["abs"])
        an.append(c["new"]["abs"])
    o = np.argsort(mus)
    mus = np.array(mus)[o]
    ao = np.array(ao)[o]
    an = np.array(an)[o]
    # extend: absorbed -> 0 at mu0 = -1 (fully dark), -> the mu0 = 1 value at mu0 = 1
    mg = np.concatenate(([-1.0], mus, [1.0]))
    og = np.concatenate(([0.0], ao, [ao[-1]/abs(mus[-1])]))
    ng = np.concatenate(([0.0], an, [an[-1]/abs(mus[-1])]))
    Io = float(np.trapezoid(og, mg))
    In = float(np.trapezoid(ng, mg))
    # pi r^2 F* at r = r(tau_slant = 1) = 1.2409e10 (measured in README.md)
    RABS = 1.2409e10
    print(f"    int(absorbed) d(mu0) over [-1,1]:  old {Io:.5e}  new {In:.5e}  "
          f"new/old {In/Io:.4f}")
    print(f"    P/(pi r_abs^2 F*) with r_abs = {RABS:.4e}: "
          f"old {2.0*np.pi*RABS**2*Io/(np.pi*RABS**2*FSTAR):.5f}  "
          f"new {2.0*np.pi*RABS**2*In/(np.pi*RABS**2*FSTAR):.5f}")
