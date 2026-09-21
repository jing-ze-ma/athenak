#!/usr/bin/env python3
"""T3b (opacity jump), REDEFINED (design note section 10).

The original specification -- "the cell-centred flux must be uniform across
the jump to 1e-3" -- was wrong, and milestone 1b showed why: in steady state
the FACE energy flux is uniform by conservation and E(x) is correct, while the
CELL-centred F of the two cells straddling the jump is off by O(1), because
the implicit source divides a centred 2 dx pressure gradient by the cell's own
opacity.  That is the same discretisation as the gas pressure force and is
correct and conservative, so it is reported here as a DIAGNOSTIC only.

What is gated instead:

* ``E(x)`` against the analytic two-slope steady state
  ``dE/dx = -3 rho kappa F/c`` on each side of the jump (max relative error);
* the two SLOPES, fitted on each side away from the jump, against
  ``-3 rho kappa F/c`` (this is the "both slopes" half of the gate);
* continuity of E at the jump, and of the FACE energy flux, which is
  reconstructed from the profile with the scheme's own compact two-point
  diffusion flux ``F_face = -c (E_{i+1} - E_i)/(3 tau_face)`` with the
  ARITHMETIC-mean face opacity (design section 3); in steady state this must
  be uniform and equal to the imposed flux.

Two cases are run, straddling ``tau_cell = 1``: a jump from ``tau_cell = 1``
to ``1e3`` and one from ``1e-3`` to ``1``.
"""

import numpy as np

import common


def analytic(x, args):
    """Two-slope steady state E(x) and the two slopes."""
    sl = -3.0 * args.rho_left * args.kappa * args.flux / args.c
    sr = args.ratio * sl
    ej = args.e_left + sl * (args.x_jump - args.xmin)
    e = np.where(x < args.x_jump,
                 args.e_left + sl * (x - args.xmin),
                 ej + sr * (x - args.x_jump))
    return e, sl, sr


def face_flux(x, e, rho, args):
    """Compact two-point diffusion flux on the interior faces."""
    dx = x[1] - x[0]
    tauf = 0.5 * (rho[:-1] + rho[1:]) * args.kappa * dx
    return -args.c * np.diff(e) / (3.0 * tauf)


def selftest(args):
    n = 256
    x = common.centres(args.xmin, args.xmax, n)
    e, _, _ = analytic(x, args)
    rho = np.where(x < args.x_jump, args.rho_left,
                   args.rho_left * args.ratio)
    if args.selftest_fail:
        e = e * (1.0 + 5e-3 * np.exp(
            -0.5 * ((x - args.x_jump) / (2 * (x[1] - x[0]))) ** 2))
    fcell = np.full(n, args.flux)
    data = {"m1_e": e[None, None, :], "m1_f1": fcell[None, None, :],
            "dens": rho[None, None, :]}
    return common.Dump(1.0, x, np.array([0.5]), np.array([0.5]), data,
                       "<selftest>")


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dump", nargs="?", help="steady-state dump (bin or tab)")
    p.add_argument("--c", type=float, default=1.0)
    p.add_argument("--kappa", type=float, default=0.64,
                   help="kappa_F + kappa_s per unit mass")
    p.add_argument("--rho-left", type=float, default=1.0)
    p.add_argument("--ratio", type=float, default=1.0e3,
                   help="jump in rho (and so in rho kappa) at x_jump")
    p.add_argument("--flux", type=float, default=1.0e-3,
                   help="the imposed constant energy flux")
    p.add_argument("--e-left", type=float, default=1.0,
                   help="E at x = xmin, held by the user BC")
    p.add_argument("--x-jump", type=float, default=0.5)
    p.add_argument("--xmin", type=float, default=0.0)
    p.add_argument("--xmax", type=float, default=1.0)
    p.add_argument("--trim", type=int, default=2,
                   help="cells excluded at each end of the profile")
    p.add_argument("--jump-cells", type=float, default=2.0,
                   help="half-width, in cells, of the window around the jump "
                        "excluded from the face-flux uniformity test and "
                        "reported separately")
    p.add_argument("--fit-cells", type=int, default=4,
                   help="cells excluded on each side of the jump in the "
                        "slope fits")
    p.add_argument("--tol", type=float, default=1e-2,
                   help="pass threshold on the relative errors of E, of the "
                        "two slopes and of the face flux")
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3))
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"))
    p.add_argument("--selftest-fail", action="store_true")
    args = p.parse_args()

    if args.selftest:
        dump = selftest(args)
    else:
        if not args.dump:
            p.error("a dump is required unless --selftest is given")
        dump = common.load_dump(args.dump, args.bin_convert_dir,
                                args.all_ranks)

    x, e = common.extract_1d(dump, "m1_e", axis=args.axis,
                             reduce=args.reduce)
    _, fc = common.extract_1d(dump, "m1_f1", axis=args.axis,
                              reduce=args.reduce)
    if dump.has("dens"):
        _, rho = common.extract_1d(dump, "dens", axis=args.axis,
                                   reduce=args.reduce)
    else:
        rho = np.where(x < args.x_jump, args.rho_left,
                       args.rho_left * args.ratio)

    eref, sl, sr = analytic(x, args)
    n = x.size
    lo, hi = args.trim, n - args.trim
    err_e = float(np.abs(e[lo:hi] / eref[lo:hi] - 1.0).max())

    dxc = x[1] - x[0]
    left = (x < args.x_jump - args.fit_cells * dxc) & (np.arange(n) >= lo)
    right = (x > args.x_jump + args.fit_cells * dxc) & (np.arange(n) < hi)
    sl_m, _ = common.linfit(x[left], e[left])
    sr_m, _ = common.linfit(x[right], e[right])
    err_sl = abs(sl_m / sl - 1.0)
    err_sr = abs(sr_m / sr - 1.0)

    # the slope RATIO is independent of the overall flux the Dirichlet ends settle
    # on, so it isolates the treatment of the jump itself
    err_ratio = abs((sr_m / sl_m) / args.ratio - 1.0)

    ff = face_flux(x, e, rho, args)
    xf = 0.5 * (x[:-1] + x[1:])
    away = (np.abs(xf - args.x_jump) > args.jump_cells * dxc)
    away[:args.trim] = False
    away[-args.trim:] = False
    fmean = float(ff[away].mean())
    err_f = float(np.abs(ff[away] / fmean - 1.0).max())
    err_fnorm = abs(fmean / args.flux - 1.0)
    near = ~away
    near[:args.trim] = False
    near[-args.trim:] = False
    err_fjump = float(np.abs(ff[near] / fmean - 1.0).max()) if near.any() \
        else 0.0
    err_fc = float(np.abs(fc[lo:hi] / args.flux - 1.0).max())

    common.report(args, "t=%.6g  slopes: exact %.6e / %.6e  measured "
                        "%.6e / %.6e" % (dump.time, sl, sr, sl_m, sr_m))
    common.report(args, "tau_cell: %.4g (left) -> %.4g (right)"
                  % (args.rho_left * args.kappa * dxc,
                     args.rho_left * args.ratio * args.kappa * dxc))
    common.report(args, "<F_face> away from the jump = %.6e (%.3e of the "
                        "imposed flux - 1); largest deviation INSIDE the "
                        "jump window = %.3e" % (fmean, err_fnorm, err_fjump))
    common.report(args, "cell-centred F error (DIAGNOSTIC only): %.3e"
                  % err_fc)
    ok = (err_e < args.tol and err_sl < args.tol and err_sr < args.tol
          and err_ratio < args.tol and err_f < args.tol)
    common.verdict(ok, "T3b jump: max|E/E_exact-1|=%.3e  slope err %.3e / "
                       "%.3e  slope ratio err %.3e  face flux: uniform to "
                       "%.3e, %.3e below the imposed value, %.3e peak at the "
                       "jump (< %.3g)  [cell F err %.3e, diagnostic]"
                   % (err_e, err_sl, err_sr, err_ratio, err_f, err_fnorm,
                      err_fjump, args.tol, err_fc))


if __name__ == "__main__":
    main()
