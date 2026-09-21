#!/usr/bin/env python3
"""T3b (opacity jump): steady constant-flux state across a 1e3 jump.

Design note section 8, (T3b) (Bloch et al. 2021 Sect. 5.2): in the steady
state the radiation flux must be uniform across a jump in rho*kappa_F,
because the face opacity is the ARITHMETIC mean of rho*kappa_F (section 3).
Pass: max |F/<F> - 1| < --tol (default 1e-3), both over the whole profile
and in the window around the jump.
"""

import numpy as np

import common


def analyse(dump, args):
    x, f = common.extract_1d(dump, args.flux_var, axis=args.axis,
                             reduce=args.reduce)
    n = x.size
    lo, hi = args.trim, n - args.trim
    xs, fs = x[lo:hi], f[lo:hi]
    fmean = float(fs.mean())
    if fmean == 0.0:
        raise ValueError("mean flux is zero")
    dev = np.abs(fs / fmean - 1.0)

    xj = args.x_jump
    if xj is None and dump.has("dens"):
        _, rho = common.extract_1d(dump, "dens", axis=args.axis,
                                   reduce=args.reduce)
        xj = float(x[int(np.argmax(np.abs(np.diff(np.log(
            np.maximum(rho, 1e-300))))))])
    dxc = x[1] - x[0]
    if xj is None:
        jdev = dev
    else:
        sel = np.abs(xs - xj) <= args.jump_cells * dxc
        jdev = dev[sel] if sel.any() else dev
    return fmean, float(dev.max()), float(jdev.max()), xj


def selftest(args):
    n = 256
    x = common.centres(0.0, 1.0, n)
    f = 1.0 + 1e-5 * np.sin(8 * np.pi * x)
    rho = np.where(x < 0.5, 1.0, 1e3)
    if args.selftest_fail:
        f = f + 5e-3 * np.exp(-0.5 * ((x - 0.5) / (2 * (x[1] - x[0]))) ** 2)
    data = {args.flux_var: f[None, None, :], "dens": rho[None, None, :]}
    return common.Dump(1.0, x, np.array([0.5]), np.array([0.5]), data,
                       "<selftest>")


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dump", nargs="?", help="steady-state bin dump")
    p.add_argument("--flux-var", default="m1_f1",
                   help="flux component along the gradient")
    p.add_argument("--x-jump", type=float, default=None,
                   help="coordinate of the opacity jump (default: from dens)")
    p.add_argument("--jump-cells", type=float, default=4.0,
                   help="half-width of the jump window, in cells")
    p.add_argument("--trim", type=int, default=2,
                   help="cells excluded at each end of the profile")
    p.add_argument("--tol", type=float, default=1e-3,
                   help="pass threshold on max |F/<F> - 1|")
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3))
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"))
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.selftest:
        dump = selftest(args)
    else:
        if not args.dump:
            p.error("a dump is required unless --selftest is given")
        dump = common.load_dump(args.dump, args.bin_convert_dir,
                                args.all_ranks)

    fmean, dmax, djump, xj = analyse(dump, args)
    common.report(args, "t=%.6g  <F>=%.6e  x_jump=%s"
                  % (dump.time, fmean, "auto-none" if xj is None
                     else "%.6g" % xj))
    ok = dmax < args.tol and djump < args.tol
    common.verdict(ok, "T3b jump: max|F/<F>-1|=%.3e overall, %.3e at the "
                       "jump (< %.3g), <F>=%.6e"
                   % (dmax, djump, args.tol, fmean))


if __name__ == "__main__":
    main()
