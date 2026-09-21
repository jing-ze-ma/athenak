#!/usr/bin/env python3
"""T4c (sheared advected pulse): convergence of the enthalpy-flux velocity.

Milestone 1c-B.  The T4 dynamic pulse is carried by a PRESCRIBED velocity
field that varies along x,

    v(x) = v0 [1 + shear_amp sin(2 pi (x - x1min)/L)],

held fixed by a user source term (see src/pgen/tests/rad_m1_tests.cpp).  The
advective split of the E equation adds back the enthalpy flux
A = v E + v.P upwinded at the face; with <rad_m1>/split_vel = cell that flux
is built from each side's CELL velocity, which is only first order at the
face once v has a gradient, and with split_vel = recon the velocity is
reconstructed there with the same method as (E, f_i).  A uniform v cannot
tell the two apart, which is why T4 could not.

There is no closed-form solution (the prescribed v is not divergence free),
so the reference is a run of the same problem on a much finer grid, block
averaged down to each coarse grid.  The script reports the relative L1 error
of E per resolution and the measured order between consecutive ones.
"""

import numpy as np

import common


def coarsen(y, factor):
    """Block average a fine profile down by an integer factor."""
    n = y.size // factor
    if n * factor != y.size:
        raise ValueError("reference length %d is not %d x the coarse grid"
                         % (y.size, factor))
    return y.reshape(n, factor).mean(axis=1)


def profile(path, all_ranks=False, bin_convert_dir=None):
    """E(x) of a 1-D m1 dump, as (x, E)."""
    dump = common.load_dump(path, bin_convert_dir=bin_convert_dir,
                            all_ranks=all_ranks)
    x, e = common.extract_1d(dump, "m1_e", axis=1)
    return np.asarray(x, dtype=float), np.asarray(e, dtype=float)


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("--ref", default=None, help="fine-grid reference dump")
    p.add_argument("--run", nargs=2, action="append", metavar=("N", "DUMP"),
                   help="resolution label and its final dump; repeat")
    p.add_argument("--max-l1", type=float, default=1.0e-2,
                   help="pass threshold on the L1 error of the finest run")
    p.add_argument("--min-order", type=float, default=0.0,
                   help="pass threshold on the measured order")
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.selftest:
        # a synthetic first-order and second-order error sequence
        ns = np.array([64.0, 128.0, 256.0])
        power = 1.0 if args.selftest_fail else 2.0
        errs = (64.0 / ns) ** power * 1.0e-3
        orders = list(np.log2(errs[:-1] / errs[1:]))
        ok = (errs[-1] <= args.max_l1) and (min(orders) >= 1.8)
        common.verdict(ok, "T4c shear(selftest): orders %s L1=%.3e"
                       % (["%.3f" % o for o in orders], errs[-1]))

    if not args.ref or not args.run:
        raise SystemExit("--ref and at least one --run are required")

    xr, er = profile(args.ref, args.all_ranks, args.bin_convert_dir)
    ns, errs = [], []
    for label, path in args.run:
        n = int(label)
        xc, ec = profile(path, args.all_ranks, args.bin_convert_dir)
        if ec.size != n:
            raise SystemExit("%s has %d cells, not %d" % (path, ec.size, n))
        eref = coarsen(er, er.size // n)
        errs.append(float(np.abs(ec - eref).sum() / np.abs(eref).sum()))
        ns.append(n)
        common.report(args, "  N=%-6d L1(E) = %.6e" % (n, errs[-1]))
    order = [float(np.log2(errs[i] / errs[i + 1]))
             for i in range(len(errs) - 1)]
    ok = (errs[-1] <= args.max_l1)
    if order:
        ok = ok and (min(order) >= args.min_order)
    common.verdict(ok, "T4c shear: L1(E) %s orders %s (N = %s)"
                   % (["%.4e" % e for e in errs],
                      ["%.3f" % o for o in order],
                      ",".join(str(n) for n in ns)))


if __name__ == "__main__":
    main()
