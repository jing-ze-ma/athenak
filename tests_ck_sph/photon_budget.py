"""Exact photon budget for the pseudo-spherical direct beam (problem/ck_beam_sph).

Gate (e) in README_beam.md reported P/(pi r_abs^2 F*) = 1.0015 (old) vs 1.1064 (new) and
called the new scheme "not photon-conserving".  This script tests that claim against the
EXACT answer for a spherically symmetric atmosphere, which is what gate (e) integrates:
one 1-D profile used for every column.  For such an atmosphere the pseudo-spherical beam
IS exact ray tracing, so the absorbed power must be

    P_exact = F* (1-albedo) * int_0^{r_top} 2 pi b a(b) db,

    a(b) = 1 - exp(-tau_chord(b))     for b > r_cut   (ray passes through, both legs)
    a(b) = 1 - exp(-tau_in(b))        for b <= r_cut  (ray dives to the opaque cut face;
                                       the remainder exp(-tau_in) is what the kernel
                                       DROPS -- see the note on the cut below)

and the scheme's global deposit is

    P_scheme = 2 pi int_{-1}^{1} dmu0 sum_i Qb_i(mu0) (r_{i+1}^3 - r_i^3)/3.

WHAT THE KERNEL DOES WITH THE POWER THAT REACHES THE CUT: it is DROPPED.
src/utils/two_stream_rt.hpp:4111 walks faces only down to icut
(`for (int i=ie; i>icut-1; --i)`) and the block ends at :4169 with no term for the
residual transmission exp(-tau(r_cut)); for a mu0 < 0 ray whose impact parameter falls
below the cut, :4121-4122 sets `dark = true` and :4150-4152 then deposits exactly zero.
Nothing is handed to the conduction/interior region.  So the correct reference for the
scheme is a(b) = 1 - exp(-tau_in(b)) for b <= r_cut, NOT a(b) = 1.  Both are reported.

PROFILES.  The correlated-k column dumps carry no opacity, so the per-g-point kappa of
gate (e) is not recoverable from the artefacts.  Two profiles are used instead:

  (A) REAL, effective grey.  The band-summed vertical transmission of the real t=0
      atmosphere is read straight out of gate (e)'s own artefact bd_p979_old/col.txt: the
      old (plane-parallel) deposit is Qb_i = (1-A) F* mu0 (T_{i+1} - T_i)/dz_i, so the
      cumulative sum of Q dz from the top gives T(r_i) exactly, hence
      tau_vert(r_i) = -mu0 ln T(r_i) and kappa rho on every shell.  This is the real
      stratification and the real total depth (T at the cut = 2.0e-4); it is grey, i.e.
      it reproduces the real vertical transmission but not the non-grey slowing of
      decay along the long twilight chords.
  (B) ANALYTIC isothermal exponential, X = r/H = 7.6 (H = 1.24e9, the production H of
      README_beam.md gate (c2)), normalised to tau_vert = 1 at r_abs = 1.2409e10; plus a
      4-g-point NON-GREY variant with the same band-mean depth, to show that the
      conservation property is per g-point and survives a non-grey mixture.

Everything (chords, weights, the tangent-shell max(), the dtau -> 0 guard, the tau > 60
break, the dark rule) is a line-by-line transcription of the kernel block at
src/utils/two_stream_rt.hpp:4102-4169; the old plane-parallel beam is the block at
:4041-4068 (tau_vert * facsw with facsw = 1/max(mu0, 0.1), dark for mu0 <= 0).
"""
import numpy as np

SIG = 5.6704e-5
TIRR = 2500.0*np.sqrt(2.0)
FSTAR = SIG*TIRR**4
ALB = 0.991541                      # (1 - albedo)
RABS = 1.2409e10                    # README.md: radius where the slant depth reaches 1
FA = FSTAR*ALB

# ----------------------------------------------------------------------- the grid ----
DUMP = "bd_p979_old/col.txt"
MU_DUMP = 9.7875341691e-01
ICUT = 19


def load_profile():
    """Faces, shell (kappa rho), pressure, over the correlated-k domain [r_cut, r_top]."""
    d = np.loadtxt(DUMP)
    idx = d[:, 0].astype(int)
    rf, p = d[:, 1], d[:, 2]
    Q = d[:, 5]
    lo = int(np.where(idx == ICUT)[0][0])
    dz = np.diff(rf)
    # cumulative absorbed above the lower face of each cell -> transmission
    S = np.concatenate((np.cumsum((Q[:-1]*dz)[::-1])[::-1], [0.0]))
    T = 1.0 - S/(FA*MU_DUMP)
    tauv = -MU_DUMP*np.log(np.clip(T, 1.0e-30, None))
    r = rf[lo:]
    kr = -np.diff(tauv[lo:])/np.diff(r)
    kr = np.maximum(kr, 0.0)
    return r, kr, p[lo:], T[lo]


def kr_exp(r, H, r_ref, tau_ref):
    """Isothermal exponential kappa rho = k0 exp(-(rc-r0)/H), tau_vert(r_ref)
    = tau_ref."""
    rc = 0.5*(r[:-1] + r[1:])
    k = np.exp(-(rc - r[0])/H)
    tv = np.cumsum((k*np.diff(r))[::-1])[::-1]          # unnormalised tau_vert at faces
    i = int(np.argmin(np.abs(r[:-1] - r_ref)))
    return k*(tau_ref/tv[i])


# ------------------------------------------------------------------ exact quadrature --
def chord_tau(b, r, kr, both_legs):
    """Optical depth of a ray of impact parameter b: full chord (both legs) if
    both_legs, else the single incoming leg from r_top down to r[0] = r_cut."""
    S = np.sqrt(np.maximum(r*r - b*b, 0.0))
    ds = np.diff(S)
    if both_legs:
        j0 = max(int(np.searchsorted(r, b, side="right") - 1), 0)
        return 2.0*float(np.sum(ds[j0:]*kr[j0:]))
    return float(np.sum(ds*kr))


def p_exact(r, kr, ngl=24):
    """P_exact / (F*(1-alb)) = int 2 pi b a(b) db, split at the cut.
    Returns (above_cut_outer, above_cut_inner, reaching_cut, opaque_body_variant)."""
    xg, wg = np.polynomial.legendre.leggauss(ngl)
    edges = np.concatenate(([0.0], r))          # [0, r_cut] then every shell face
    inner = outer = drop = 0.0
    for a, c in zip(edges[:-1], edges[1:]):
        bb = 0.5*(c - a)*xg + 0.5*(c + a)
        ww = 0.5*(c - a)*wg
        for b, w in zip(bb, ww):
            if b <= r[0]:
                t = chord_tau(b, r, kr, False)
                inner += w*2.0*np.pi*b*(1.0 - np.exp(-t))
                drop += w*2.0*np.pi*b*np.exp(-t)
            else:
                t = chord_tau(b, r, kr, True)
                outer += w*2.0*np.pi*b*(1.0 - np.exp(-t))
    return outer, inner, drop, outer + inner + drop


# ------------------------------------------------------------------- the two kernels --
def qb_new(mu0, r, krs, wts=None):
    """src/utils/two_stream_rt.hpp:4102-4169, transcribed.  krs: (ng, nshell)."""
    n = r.size - 1
    ng = krs.shape[0]
    if wts is None:
        wts = np.ones(ng)
    sinz = np.sqrt(max(0.0, 1.0 - mu0*mu0))
    rcut = r[0]
    Q = np.zeros(n)
    tauh = np.zeros(ng)
    thi = np.ones(ng)
    for i in range(n - 1, -1, -1):
        b = r[i]*sinz
        b2 = b*b
        dark = (mu0 < 0.0) and (b <= rcut)
        if not dark:
            S = np.sqrt(np.maximum(r*r - b2, 0.0))
            ds = np.diff(S)
            w = np.ones(n)
            if mu0 < 0.0:
                jlo = max(int(np.searchsorted(r, b, side="right") - 1), 0)
                w[:i] = 2.0
            else:
                jlo = i
            taul = (ds[jlo:]*w[jlo:]*krs[:, jlo:]).sum(axis=1)
        else:
            taul = np.full(ng, 1.0e30)
        dtl = taul - tauh
        tlo = np.where(dark, 0.0, np.exp(-np.minimum(taul, 700.0)))
        fac = np.where(dtl > 1.0e-3, (thi - tlo)/np.where(dtl > 1.0e-3, dtl, 1.0),
                       thi*(1.0 - 0.5*dtl))
        Q[i] = float(np.sum(wts*fac*krs[:, i]))
        tauh, thi = taul, tlo
        if tauh.min() > 60.0:
            break
    return Q


def qb_old(mu0, r, krs, wts=None):
    """The plane-parallel beam, src/utils/two_stream_rt.hpp:4041-4068."""
    n = r.size - 1
    if wts is None:
        wts = np.ones(krs.shape[0])
    Q = np.zeros(n)
    if mu0 <= 0.0:
        return Q
    facsw = 1.0/mu0 if mu0 > 0.1 else 10.0
    dz = np.diff(r)
    tau = np.zeros(krs.shape[0])
    tr = np.ones(krs.shape[0])
    for i in range(n - 1, -1, -1):
        tau = tau + krs[:, i]*dz[i]
        tn = np.exp(-np.minimum(tau*facsw, 700.0))
        Q[i] = float(np.sum(wts*(tr - tn)))/facsw/dz[i]
        tr = tn
    return Q


def p_scheme(r, krs, kern, nmu=400, wts=None):
    """2 pi int dmu0 sum_i Q_i dV_i / (F*(1-alb)), dV per steradian."""
    xg, wg = np.polynomial.legendre.leggauss(nmu)
    dv = (r[1:]**3 - r[:-1]**3)/3.0
    tot = 0.0
    for mu0, w in zip(xg, wg):
        tot += w*float(np.sum(kern(mu0, r, krs, wts)*dv))
    return 2.0*np.pi*tot


# -------------------------------------------------------------------------- report ----
def budget(tag, r, krs, note="", wts=None):
    print(f"\n=== {tag} ===  {note}")
    ng = krs.shape[0]
    if wts is None:
        wts = np.ones(ng)
    if ng == 1:
        out, inn, drop, _ = p_exact(r, krs[0])
    else:                       # non-grey: the budget is linear in the g-points
        out = inn = drop = 0.0
        for g in range(ng):
            o, i2, d, _ = p_exact(r, krs[g])
            out += wts[g]*o
            inn += wts[g]*i2
            drop += wts[g]*d
    Pex = out + inn
    Pfull = out + inn + drop        # if the body below the cut absorbed the residual
    Pn = p_scheme(r, krs, qb_new, wts=wts)
    Po = p_scheme(r, krs, qb_old, wts=wts)
    reff = np.sqrt(Pex/np.pi)
    print(f"  P_exact (scheme-consistent: residual at the cut dropped) "
          f"= {Pex:.6e} cm^2 x F*(1-alb)")
    print(f"    of which absorbed at b >  r_cut (limb/twilight chords) = {out:.6e} "
          f"({100*out/Pex:.2f} %)")
    print(f"    of which absorbed at b <= r_cut (dayside, above the cut) = {inn:.6e} "
          f"({100*inn/Pex:.2f} %)")
    print(f"    power REACHING the cut face and dropped by the kernel  = {drop:.6e} "
          f"({100*drop/Pex:.3f} % of P_exact)")
    print(f"    (if the body absorbed it instead: P = {Pfull:.6e}, "
          f"r_eff = {np.sqrt(Pfull/np.pi):.5e})")
    print(f"  P_new    = {Pn:.6e}   P_new/P_exact = {Pn/Pex:.6f}")
    print(f"  P_old    = {Po:.6e}   P_old/P_exact = {Po/Pex:.6f}")
    print(f"  r_eff = sqrt(P_exact/pi) = {reff:.5e} cm    "
          f"P_exact/(pi r_abs^2) = {Pex/(np.pi*RABS**2):.5f}")
    return Pex, Pn, Po, reff


def radii(r, kr, p=None):
    """r(tau_vert = 1) and the impact parameter with full-chord tau = 1."""
    tv = np.concatenate((np.cumsum((kr*np.diff(r))[::-1])[::-1], [0.0]))
    r1 = float(np.interp(0.0, np.log(np.clip(tv[::-1], 1e-300, None)), r[::-1]))
    bs = np.linspace(r[0], r[-1], 20001)
    tc = np.array([chord_tau(b, r, kr, True) for b in bs])
    rl = float(np.interp(0.0, np.log(np.clip(tc[::-1], 1e-300, None)), bs[::-1]))
    return r1, rl


if __name__ == "__main__":
    r, kr, p, Tcut = load_profile()
    print("EXACT PHOTON BUDGET for problem/ck_beam_sph")
    print(f"  grid: {r.size-1} shells, r_cut = {r[0]:.5e}, r_top = {r[-1]:.5e} cm")
    print(f"  F* = {FSTAR:.5e} erg/s/cm2, (1-albedo) = {ALB}")
    print(f"  (A) real effective-grey profile from {DUMP}: band-summed vertical "
          f"transmission at the cut = {Tcut:.4e}, tau_vert(cut) = "
          f"{-MU_DUMP*np.log(Tcut):.4f}")

    # pressure scale height near the tau = 1 level
    rc = 0.5*(r[:-1] + r[1:])
    H = -np.diff(r)/np.log(p[1:]/p[:-1])
    ia = int(np.argmin(np.abs(rc - RABS)))
    Hab = float(H[ia])
    print(f"  pressure scale height at r_abs: H = {Hab:.4e} cm  (X = r/H = "
          f"{RABS/Hab:.2f})")

    PA = budget("(A) REAL effective-grey profile", r, kr[None, :],
                "opacity reconstructed from the old scheme's own vertical transmission")
    r1, rl = radii(r, kr)
    print(f"  r(tau_vert = 1)      = {r1:.5e} cm")
    print(f"  r(tau_chord,limb = 1) = {rl:.5e} cm   "
          f"= r(tau_vert=1) + {(rl-r1)/Hab:.2f} H")
    print(f"  r_eff                = {PA[3]:.5e} cm   "
          f"= r(tau_vert=1) + {(PA[3]-r1)/Hab:.2f} H")
    print(f"  r_abs (README)       = {RABS:.5e} cm")
    print(f"  (pi r_eff^2)/(pi r_abs^2) = {(PA[3]/RABS)**2:.5f}")

    # (B) analytic
    H0 = 1.24e9
    k1 = kr_exp(r, H0, RABS, 1.0)
    budget("(B1) ANALYTIC isothermal exponential, X = 7.6", r, k1[None, :],
           f"H = {H0:.3e}, tau_vert(r_abs) = 1")
    # non-grey: 4 g-points spanning 3 decades in kappa, equal weights
    fac = np.array([0.03, 0.3, 3.0, 30.0])
    kng = np.array([f*k1 for f in fac])
    budget("(B2) ANALYTIC, NON-GREY 4 g-points (kappa x 0.03 .. 30)", r, kng,
           "equal weights 1/4", wts=np.full(4, 0.25))

    # mu0-quadrature convergence of the scheme integral
    print("\n=== mu0-quadrature convergence of P_new, profile (A) ===")
    for nmu in [50, 100, 200, 400, 800]:
        print(f"  nmu = {nmu:4d}: P_new = {p_scheme(r, kr[None, :], qb_new, nmu):.6e}")

    # ---- validation of the transcription against the code's own column dumps, and the
    # decomposition of gate (e)'s 1.1064
    print("\n=== transcription vs the code, and gate (e)'s quadrature ===")
    print("  'absorbed per unit top area' = sum_i Q_i dz_i, as in beamprof.py")
    print("    mu0    |   old: code / this script   |   new: code / this script")
    tags = [("p979", 9.7875341691e-01), ("p435", 4.3472405714e-01),
            ("p102", 1.0198e-01), ("n102", -1.0198e-01),
            ("n195", -1.9470e-01), ("n513", -5.1280e-01)]
    dzc = np.diff(r)
    mus, co, cn, so, sn = [], [], [], [], []
    for tag, mu in tags:
        mu_r = mu
        vals = {}
        for w, kern in (("old", qb_old), ("new", qb_new)):
            d = np.loadtxt(f"bd_{tag}_{w}/col.txt")
            with open(f"bd_{tag}_{w}/col.txt") as f:
                for ln in f:
                    if "mu0 =" in ln:
                        mu_r = float(ln.split("mu0 =")[1].split(",")[0])
                        break
            lo = int(np.where(d[:, 0].astype(int) == ICUT)[0][0])
            vals[w] = (float(np.sum(d[lo:-1, 5]*np.diff(d[lo:, 1]))),
                       FA*float(np.sum(kern(mu_r, r, kr[None, :])*dzc)))
        mus.append(mu_r)
        co.append(vals["old"][0])
        cn.append(vals["new"][0])
        so.append(vals["old"][1])
        sn.append(vals["new"][1])
        print(f"  {mu_r:+7.4f} | {vals['old'][0]:11.4e} / {vals['old'][1]:11.4e}"
              f"  ({vals['old'][1]/max(vals['old'][0], 1e-30):6.3f}) |"
              f" {vals['new'][0]:11.4e} / {vals['new'][1]:11.4e}"
              f"  ({vals['new'][1]/max(vals['new'][0], 1e-30):6.3f})")
    o = np.argsort(mus)
    mg = np.concatenate(([-1.0], np.array(mus)[o], [1.0]))
    for lab, a in (("old (code)", np.array(co)[o]), ("new (code)", np.array(cn)[o])):
        g = np.concatenate(([0.0], a, [a[-1]/abs(np.array(mus)[o][-1])]))
        Ig = float(np.trapezoid(g, mg))
        print(f"  gate (e) recipe, {lab}: P/(pi r_abs^2 F*) = "
              f"{2.0*Ig/FSTAR:.5f}")
    print("  (gate (e) puts every photon at r_abs and samples mu0 at 6 points; the"
          " correct\n   weight is r^2 dr inside the column, which is what P_new above"
          " uses.)")

    # the real (non-grey) atmosphere with the CORRECT radial weight, still only the six
    # sampled mu0: sum_i Q_i (r_{i+1}^3 - r_i^3)/3 instead of sum_i Q_i dz_i r_abs^2
    print("\n=== the code's own columns with the correct r^2 dr weight ===")
    vo, vn = [], []
    for tag, _ in tags:
        for w, acc in (("old", vo), ("new", vn)):
            d = np.loadtxt(f"bd_{tag}_{w}/col.txt")
            lo = int(np.where(d[:, 0].astype(int) == ICUT)[0][0])
            rr = d[lo:, 1]
            acc.append(float(np.sum(d[lo:-1, 5]*(rr[1:]**3 - rr[:-1]**3)/3.0)))
    for lab, a in (("old", np.array(vo)[o]), ("new", np.array(vn)[o])):
        g = np.concatenate(([0.0], a, [a[-1]/abs(np.array(mus)[o][-1])]))
        Ig = float(np.trapezoid(g, mg))
        print(f"  {lab}: P/(pi r_abs^2 F*) = {2.0*Ig/(RABS**2*FSTAR):.5f}"
              f"   r_eff = {np.sqrt(2.0*Ig/FA):.5e} cm")
