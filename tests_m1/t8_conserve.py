#!/usr/bin/env python3
"""T8c: conservation of the coupled gas + radiation totals.

Design note sections 1 and 8.  With ``chat`` on the radiation time
derivative only, what the scheme conserves over a periodic box is

    E_tot = sum [ e_int + (1/2) rho v^2 + (c/chat) E ]
    P_tot = sum [ rho v_x + F_x/(chat c) ]

(at chat = c these are the true totals).  The script reads a series of
``file_type = tab`` dumps -- the bin writer is single precision and cannot
resolve the drifts this gate targets -- pairs each ``*.m1.*`` dump with the
``*.hydro_w.*`` dump of the same index, and reports the relative drift of
both totals from the first dump to the last.

``P_tot`` is reported as an ABSOLUTE drift normalised by sum|rho v_x| +
sum|F_x/(chat c)|, because the momentum of a symmetric pulse in a uniformly
advected box is not near zero but its two contributions can cancel.
"""

import numpy as np

import common


def totals(m1_path, hyd_path, c, chat, all_ranks=False, bin_convert_dir=None):
    dm = common.load_dump(m1_path, bin_convert_dir=bin_convert_dir,
                          all_ranks=all_ranks)
    dh = common.load_dump(hyd_path, bin_convert_dir=bin_convert_dir,
                          all_ranks=all_ranks)
    _, e = common.extract_1d(dm, "m1_e", axis=1)
    _, f1 = common.extract_1d(dm, "m1_f1", axis=1)
    _, d = common.extract_1d(dh, "dens", axis=1)
    _, vx = common.extract_1d(dh, "velx", axis=1)
    _, ei = common.extract_1d(dh, "eint", axis=1)
    e = np.asarray(e, float)
    f1 = np.asarray(f1, float)
    d = np.asarray(d, float)
    vx = np.asarray(vx, float)
    ei = np.asarray(ei, float)
    egas = ei + 0.5 * d * vx * vx
    etot = float(egas.sum() + (c / chat) * e.sum())
    pgas = d * vx
    prad = f1 / (chat * c)
    ptot = float(pgas.sum() + prad.sum())
    pscale = float(np.abs(pgas).sum() + np.abs(prad).sum())
    return dm.time, etot, ptot, pscale


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("m1", nargs="*", help="the *.m1.*.tab dumps, in order")
    p.add_argument("--hydro", nargs="*", default=None,
                   help="the matching *.hydro_w.*.tab dumps")
    p.add_argument("--c", type=float, default=2.99792458e10)
    p.add_argument("--chat", type=float, default=None,
                   help="reduced light speed (default: --c)")
    p.add_argument("--tol", type=float, default=1.0e-12,
                   help="pass threshold on both relative drifts")
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.selftest:
        de = 3.0e-16 if not args.selftest_fail else 1.0e-6
        ok = (abs(de) <= args.tol)
        common.verdict(ok, "T8c conserve(selftest): dE/E=%.3e dP=%.3e"
                       % (de, de))

    chat = args.chat if args.chat else args.c
    if not args.m1 or not args.hydro:
        raise SystemExit("need the m1 dumps and --hydro dumps")
    if len(args.m1) != len(args.hydro):
        raise SystemExit("m1 and --hydro lists differ in length")

    rows = [totals(a, b, args.c, chat, args.all_ranks, args.bin_convert_dir)
            for a, b in zip(sorted(args.m1), sorted(args.hydro))]
    t0, e0, p0, ps0 = rows[0]
    de_max, dp_max = 0.0, 0.0
    for t, e, pp, ps in rows:
        de = (e - e0) / abs(e0)
        dp = (pp - p0) / max(ps0, 1.0e-300)
        de_max = max(de_max, abs(de))
        dp_max = max(dp_max, abs(dp))
        common.report(args, "  t=%.6e  dE/E=%+.6e  dP/|P|=%+.6e" % (t, de, dp))
    ok = (de_max <= args.tol) and (dp_max <= args.tol)
    common.verdict(ok, "T8c conserve: max|dE/E|=%.3e  max|dP|/|P|=%.3e "
                       "(<= %.1e, %d dumps)"
                   % (de_max, dp_max, args.tol, len(rows)))


if __name__ == "__main__":
    main()
