#!/usr/bin/env python3
"""T1 (free streaming): convergence of a 1-D translating pulse.

kappa = 0, so the exact solution of dE/dt + c dE/dx = 0 is the initial
profile translated by c*(t - t0) (periodic box).  For each resolution the
script computes the relative L1 error of the final dump against the
translated INITIAL dump, then the measured order between consecutive
resolutions.  Pass: order >= --min-order (default 1.8, i.e. PLM second
order; design note section 3).
"""

import numpy as np

import common


def periodic_shift(x, y, shift, xl, xr):
    """Evaluate the periodic profile y(x) at x - shift."""
    length = xr - xl
    xe = np.concatenate([x - length, x, x + length])
    ye = np.concatenate([y, y, y])
    xq = (x - shift - xl) % length + xl
    return np.interp(xq, xe, ye)


def one_resolution(init, final, args):
    x, e0 = common.extract_1d(init, "m1_e", axis=args.axis,
                              reduce=args.reduce)
    _, e1 = common.extract_1d(final, "m1_e", axis=args.axis,
                              reduce=args.reduce)
    dx = x[1] - x[0]
    xl = x[0] - 0.5 * dx
    xr = x[-1] + 0.5 * dx
    shift = args.c * (final.time - init.time)
    exact = periodic_shift(x, e0, shift, xl, xr)
    return x.size, common.l1_rel(e1, exact), shift


def selftest(args):
    """Synthetic: error amplitude ~ N^-2 (pass) or ~ N^-1 (fail)."""
    res = []
    order = 1.0 if args.selftest_fail else 2.0
    for n in (64, 128, 256):
        x = common.centres(0.0, 1.0, n)
        e0 = np.exp(-0.5 * ((x - 0.3) / 0.05) ** 2)
        shift = args.c * 0.1
        exact = periodic_shift(x, e0, shift, 0.0, 1.0)
        err = 0.5 * (64.0 / n) ** order * np.sin(2 * np.pi * x)
        res.append((n, common.l1_rel(exact + err, exact), shift))
    return res


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("--run", nargs=3, action="append", metavar=("N", "INIT",
                                                               "FINAL"),
                   help="resolution label and the initial/final dumps; "
                        "repeat for every resolution")
    p.add_argument("--c", type=float, default=1.0,
                   help="speed of light in code units (match the athinput)")
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3),
                   help="propagation axis")
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"),
                   help="how to collapse the transverse directions")
    p.add_argument("--min-order", type=float, default=1.8,
                   help="pass threshold on the measured order")
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.selftest:
        res = selftest(args)
    else:
        if not args.run or len(args.run) < 2:
            p.error("at least two --run N INIT FINAL triples are needed")
        res = []
        for _, ipath, fpath in args.run:
            init = common.load_dump(ipath, args.bin_convert_dir,
                                    args.all_ranks)
            final = common.load_dump(fpath, args.bin_convert_dir,
                                     args.all_ranks)
            res.append(one_resolution(init, final, args))
    res.sort(key=lambda r: r[0])

    for n, err, shift in res:
        common.report(args, "N=%5d  L1=%.6e  shift=%.6g" % (n, err, shift))
    orders = []
    for (n0, e0, _), (n1, e1, _) in zip(res[:-1], res[1:]):
        if e1 <= 0.0 or e0 <= 0.0:
            orders.append(float("inf"))
        else:
            orders.append(np.log(e0 / e1) / np.log(float(n1) / float(n0)))
        common.report(args, "order(%d->%d) = %.3f" % (n0, n1, orders[-1]))
    omin = min(orders)
    common.verdict(omin >= args.min_order,
                   "T1 pulse: orders [%s], min=%.3f (>= %.3g), "
                   "L1(finest)=%.3e"
                   % (", ".join("%.3f" % o for o in orders), omin,
                      args.min_order, res[-1][1]))


if __name__ == "__main__":
    main()
