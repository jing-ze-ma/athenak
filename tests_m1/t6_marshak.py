#!/usr/bin/env python3
"""T6 (Marshak wave / Su-Olson): non-equilibrium benchmark.

PROVENANCE, read this before quoting a number.  The published Su & Olson
(1996, JQSRT 56, 337) tables could not be accessed from this machine, and
their semi-analytic quadrature could not be transcribed from the paper
either.  No table values are invented here.  Instead this script builds its
OWN reference by solving the same linearised problem with a converged
discrete-ordinates (S_N) transport solver:

    (1/c) dI/dt + mu dI/dx + sigma I = (sigma c V + Q) / 2
    dV/dt = eps sigma c (E - V),   E = (1/c) int I dmu,  V = a T^4

which is the Su-Olson system: grey, isotropic, constant opacity, a material
heat capacity c_v = alpha T^3 (so that V = eps * U_material, with
eps = 4 a / alpha), a unit source in 0 < x < x0 switched off at t = t0, and a
cold medium initially.  The solver is second order in space (minmod-limited
upwind per ordinate), SSP-RK2 in time, Gauss-Legendre in mu; --ref-
convergence refines nx and nmu and reports the self-convergence, which is
the only accuracy claim made.

If you DO have the published table, feed it with --table FILE (whitespace
columns: x  E  V, already converted to code units) and the solver is
bypassed.  Pass: L1 error of E and of the material energy <= --tol
(default 2 %, design note section 8 T6: Bloch et al.'s AP scheme reaches
1.1 %, uncorrected HLL 84 %).
"""

import numpy as np

import common


def minmod(a, b):
    return np.where(a * b > 0.0, np.where(np.abs(a) < np.abs(b), a, b), 0.0)


def sn_reference(args, nx=None, nmu=None):
    """Converged S_N solution, in code units.

    Two problems, selected by ``--marshak``:

    * default: the SU-OLSON system (``c_v = alpha T^3``, so ``V = a T^4`` is
      linear in the material energy), a volumetric source in ``0 < x < x0``
      switched off at ``t0``, reflecting at ``x = 0``;
    * ``--marshak``: the CONSTANT-``c_v`` non-equilibrium Marshak wave, which
      is what an ideal-gas EOS can represent (design note section 10: ``e``
      proportional to ``T^4`` would need a user EOS).  The material energy
      ``U = rho c_v T`` is evolved with ``dU/dt = sigma c (E - a T^4)``, the
      emission term is nonlinear, and the left boundary is an incident
      isotropic bath of energy density ``--e-bath`` (``I(mu>0) = c E_b/2``,
      i.e. an incoming flux ``c E_b/4``) instead of a reflection.

    Returns (x_centres, E(x), V(x)) at t = --tend on [0, Lz]; ``V`` is
    ``a T^4`` in both cases, which is what the code's dumps are converted to.
    """
    nx = nx or args.nx
    nmu = nmu or args.nmu
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    x = common.centres(0.0, args.lz, nx)
    dx = args.lz / nx
    sig = args.sigma
    cee = args.c
    mk = args.marshak
    inten = np.zeros((nmu, nx))
    if mk:
        umat = np.full(nx, args.rho * args.cv * args.t_init)
        vmat = args.a_rad * args.t_init ** 4 * np.ones(nx)
        inten = 0.5 * cee * vmat[None, :] * np.ones((nmu, 1))
        iin = 0.5 * cee * args.e_bath
    else:
        vmat = np.zeros(nx)
        umat = None
    qsrc = np.where(x < args.x0, args.q0, 0.0)
    dt = args.cfl * dx / cee
    nt = max(1, int(np.ceil(args.tend / dt)))
    dt = args.tend / nt

    def rhs(inten, mat, time):
        big = np.zeros((nmu, nx + 4))
        big[:, 2:-2] = inten
        if mk:
            big[:, 1] = iin
            big[:, 0] = iin
            big[:, -2] = inten[:, -1]
            big[:, -1] = inten[:, -1]
        else:
            big[:, 1] = inten[::-1, 0]
            big[:, 0] = inten[::-1, 1]
        d1 = big[:, 1:-1] - big[:, :-2]
        d2 = big[:, 2:] - big[:, 1:-1]
        slope = minmod(d1, d2)
        lstate = big[:, 1:-2] + 0.5 * slope[:, 0:-1]
        rstate = big[:, 2:-1] - 0.5 * slope[:, 1:]
        up = np.where(mu[:, None] > 0.0, lstate, rstate)
        flux = cee * mu[:, None] * up
        didt = -(flux[:, 1:] - flux[:, :-1]) / dx
        erad = (wmu[:, None] * inten).sum(axis=0) / cee
        if mk:
            tt = mat / (args.rho * args.cv)
            vv = args.a_rad * tt ** 4
            didt += cee * sig * (0.5 * cee * vv[None, :] - inten)
            dmdt = sig * cee * (erad - vv)
        else:
            qnow = qsrc if time < args.t0 else 0.0 * qsrc
            didt += cee * sig * (0.5 * cee * mat[None, :] - inten)
            didt += 0.5 * cee * qnow[None, :]
            dmdt = args.eps * sig * cee * (erad - mat)
        return didt, dmdt

    mat = umat if mk else vmat
    time = 0.0
    for _ in range(nt):
        k1i, k1v = rhs(inten, mat, time)
        i1 = inten + dt * k1i
        v1 = mat + dt * k1v
        k2i, k2v = rhs(i1, v1, time + dt)
        inten = 0.5 * (inten + i1 + dt * k2i)
        mat = 0.5 * (mat + v1 + dt * k2v)
        time += dt
    erad = (wmu[:, None] * inten).sum(axis=0) / cee
    if mk:
        vmat = args.a_rad * (mat / (args.rho * args.cv)) ** 4
    else:
        vmat = mat
    return x, erad, vmat


def load_table(path):
    dat = np.loadtxt(path)
    return dat[:, 0], dat[:, 1], dat[:, 2]


def code_profiles(dump, args):
    x, erad = common.extract_1d(dump, "m1_e", axis=args.axis,
                                reduce=args.reduce)
    xabs = np.abs(x - args.centre)
    if args.material_var and dump.has(args.material_var):
        _, umat = common.extract_1d(dump, args.material_var, axis=args.axis,
                                    reduce=args.reduce)
        if args.material_is_v:
            vmat = umat
        elif args.marshak:
            # constant c_v: the dump holds the internal energy density, and
            # the reference is compared in a T^4
            if dump.has("dens"):
                _, rho = common.extract_1d(dump, "dens", axis=args.axis,
                                           reduce=args.reduce)
            else:
                rho = args.rho
            vmat = args.a_rad * (umat / (rho * args.cv)) ** 4
        else:
            vmat = args.eps * umat
    else:
        vmat = None
    return xabs, erad, vmat


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dump", nargs="?", help="dump at t = --tend (bin or tab)")
    p.add_argument("--hydro", default=None,
                   help="the matching hydro dump (file_type = tab writes one "
                        "file per <output> block)")
    p.add_argument("--table", default=None,
                   help="published reference table (x E V, code units); "
                        "bypasses the built-in S_N solver")
    p.add_argument("--sigma", type=float, default=1.0,
                   help="rho * kappa, per unit length")
    p.add_argument("--eps", type=float, default=1.0,
                   help="eps = 4 a / alpha (Su-Olson retardation parameter)")
    p.add_argument("--marshak", action="store_true",
                   help="constant-c_v non-equilibrium Marshak wave with an "
                        "incident bath at x = 0, instead of Su-Olson")
    p.add_argument("--cv", type=float, default=1.5,
                   help="--marshak: material heat capacity per unit mass "
                        "(an ideal gas has c_v = 1/(gamma-1) in code units)")
    p.add_argument("--rho", type=float, default=1.0,
                   help="--marshak: material density")
    p.add_argument("--a-rad", type=float, default=1.0,
                   help="--marshak: radiation constant in code units")
    p.add_argument("--t-init", type=float, default=0.0,
                   help="--marshak: initial material/radiation temperature")
    p.add_argument("--e-bath", type=float, default=1.0,
                   help="--marshak: energy density of the incident bath")
    p.add_argument("--c", type=float, default=1.0,
                   help="speed of light in code units")
    p.add_argument("--x0", type=float, default=0.5,
                   help="half width of the source region")
    p.add_argument("--t0", type=float, default=1e30,
                   help="time at which the source is switched off")
    p.add_argument("--q0", type=float, default=1.0,
                   help="source strength (energy density per unit time)")
    p.add_argument("--tend", type=float, default=1.0,
                   help="time of the dump (must match it)")
    p.add_argument("--lz", type=float, default=10.0,
                   help="length of the reference half-domain")
    p.add_argument("--nx", type=int, default=2000,
                   help="reference cells on [0, lz]")
    p.add_argument("--nmu", type=int, default=16,
                   help="Gauss-Legendre ordinates")
    p.add_argument("--cfl", type=float, default=0.4)
    p.add_argument("--centre", type=float, default=0.0,
                   help="symmetry centre of the run in x1")
    p.add_argument("--xmax", type=float, default=None,
                   help="compare only |x - centre| <= this")
    p.add_argument("--material-var", default="eint",
                   help="dump variable holding the material energy density")
    p.add_argument("--material-is-v", action="store_true",
                   help="that variable is already a T^4, not U_material")
    p.add_argument("--tol", type=float, default=0.02,
                   help="pass threshold on the relative L1 errors")
    p.add_argument("--ref-convergence", action="store_true",
                   help="refine the reference and report its self-error")
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3))
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"))
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.table:
        xr, er, vr = load_table(args.table)
        src = "table %s" % args.table
    else:
        xr, er, vr = sn_reference(args)
        src = "built-in S_N (nx=%d, nmu=%d)" % (args.nx, args.nmu)

    if args.ref_convergence and not args.table:
        x2, e2, v2 = sn_reference(args, nx=2 * args.nx, nmu=2 * args.nmu)
        ei = np.interp(xr, x2, e2)
        vi = np.interp(xr, x2, v2)
        common.report(args, "reference self-error: E %.3e  V %.3e"
                      % (common.l1_rel(er, ei), common.l1_rel(vr, vi)))

    if args.selftest:
        xc = xr.copy()
        ec = er.copy()
        vc = vr.copy()
        wig = 1.0 + 0.005 * np.sin(6.0 * np.pi * xc / max(xc[-1], 1e-30))
        ec = ec * wig
        vc = vc * wig
        if args.selftest_fail:
            ec = ec * 1.05
            vc = vc * 1.05
        time = args.tend
    else:
        if not args.dump:
            p.error("a dump is required unless --selftest is given")
        dump = common.load_dump(args.dump, args.bin_convert_dir,
                                args.all_ranks)
        if args.hydro:
            hyd = common.load_dump(args.hydro, args.bin_convert_dir,
                                   args.all_ranks)
            for k, v in hyd.data.items():
                dump.data.setdefault(k, v)
        xc, ec, vc = code_profiles(dump, args)
        time = dump.time

    sel = np.ones_like(xc, dtype=bool)
    if args.xmax is not None:
        sel = xc <= args.xmax
    eref = np.interp(xc[sel], xr, er)
    err_e = common.l1_rel(ec[sel], eref)
    err_v = float("nan")
    if vc is not None:
        vref = np.interp(xc[sel], xr, vr)
        err_v = common.l1_rel(np.asarray(vc)[sel], vref)

    common.report(args, "reference: %s  t=%.6g  npts=%d"
                  % (src, time, int(sel.sum())))
    ok = err_e <= args.tol and (np.isnan(err_v) or err_v <= args.tol)
    common.verdict(ok, "T6 marshak: L1(E)=%.4f  L1(material)=%.4f "
                       "(<= %.3g)  reference=%s"
                   % (err_e, err_v, args.tol, src))


if __name__ == "__main__":
    main()
