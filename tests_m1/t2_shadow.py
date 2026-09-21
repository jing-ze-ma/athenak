#!/usr/bin/env python3
"""T2, the Hayes & Norman / HERACLES shadow test (design note section 8).

A dense elliptical clump with a smooth edge sits in a thin ambient medium
and is illuminated from the inner-x1 face by a source at ``T_r``.  Two
numbers are gated on the time series of 2-D dumps:

* the ARRIVAL TIME of the radiation front at a probe far downstream of the
  clump, in the lit lane, against ``x_probe/c`` (the front is free
  streaming there, ``tau`` across the ambient medium being 0.1).  It is
  found by linear interpolation in time of ``E(x_probe, y_lit)`` through
  ``--arrival-frac`` of the incident energy density.  The default
  threshold is 1 %%, i.e. the LEADING EDGE: the half-maximum of a front
  that PLM + HLL have smeared over tens of cells lags by about half that
  width, which is a reconstruction property and not a wave speed.  The
  fitted front speed is reported at whatever threshold is chosen;
* the SHADOW DEPTH, ``E(x_shadow, y_shadow)/E(x_shadow, y_lit)`` at the
  final time, with the whole time series reported as well.  M1 shadows are
  known to stay sharp; a diffusion closure fills them in, which is what
  ``--control`` (the final dump of an otherwise identical
  ``closure = eddington`` run) exhibits.

The ambient medium is not perfectly transparent, so the lit value is read
from the run rather than assumed, and the depth is a RATIO of two cells at
the same x.
"""

import numpy as np

import common


def probe(dump, x, y):
    """E at (x, y) of a 2-D dump, bilinear (common.sample_2d)."""
    return float(common.sample_2d(dump.x1, dump.x2,
                                  np.asarray(dump.var("m1_e"))[0], x, y)[0])


def ycut(dump, x):
    """E(y) at fixed x: the nearest column of the dump."""
    i = int(np.argmin(np.abs(dump.x1 - x)))
    arr = dump.data["m1_e"]
    return np.asarray(dump.x2), np.asarray(arr[0, :, i])


def selftest(args):
    """A free-streaming front plus a perfectly sharp (or, for the failing
    variant, a filled-in) geometric shadow, on the real grid."""
    nx, ny = 140, 40
    x = common.centres(0.0, 1.0, nx)
    y = common.centres(0.0, 0.12, ny)
    xx, yy = np.meshgrid(x, y)
    dumps = []
    e_amb = args.arad * args.t0 ** 4
    for it in range(1, 41):
        t = it * args.t_cross / 20.0
        lit = (xx < args.c * t)
        # geometric shadow of the ellipse for a plane-parallel beam in +x
        shade = (yy < args.shadow_ay) & (xx > args.shadow_x0)
        e = np.where(lit, args.e_in, e_amb)
        if args.selftest_fail:
            # a diffusion closure fills the shadow in
            e = np.where(shade & lit, 0.5 * args.e_in, e)
        else:
            e = np.where(shade & lit, e_amb, e)
        data = {"m1_e": e[None, :, :]}
        dumps.append(common.Dump(t, x, y, np.array([0.5]), data,
                                 "<selftest%d>" % it))
    return dumps


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dumps", nargs="*", help="the 2-D time series (bin)")
    p.add_argument("--c", type=float, default=2.99792458e10)
    p.add_argument("--arad", type=float, default=7.5657e-15)
    p.add_argument("--t0", type=float, default=290.0,
                   help="ambient (and gas) temperature")
    p.add_argument("--tr", type=float, default=1740.0,
                   help="temperature of the illuminating source")
    p.add_argument("--x-arrival", type=float, default=0.98,
                   help="x of the arrival-time probe")
    p.add_argument("--x-shadow", type=float, default=0.8,
                   help="x of the shadow-depth probe")
    p.add_argument("--y-shadow", type=float, default=0.0)
    p.add_argument("--y-lit", type=float, default=0.115,
                   help="y of the LIT reference lane.  The design note names "
                        "0.1, but the Fermi edge of the clump (delta = 10) "
                        "puts rho = 7.7 rho0 at y = 0.09 and 2.3 rho0 at "
                        "y = 0.1, so y = 0.1 sits in the halo's penumbra; "
                        "0.115 is the first fully lit lane")
    p.add_argument("--shadow-x0", type=float, default=0.5)
    p.add_argument("--shadow-ay", type=float, default=0.06)
    p.add_argument("--arrival-frac", type=float, default=0.01,
                   help="threshold, as a fraction of the incident E, that "
                        "defines the front.  The default 0.01 is the LEADING "
                        "EDGE, which is what travels at c; at 0.5 the "
                        "PLM+HLL smearing of the front adds ~8 %% of lag")
    p.add_argument("--t-cross", type=float, default=None,
                   help="expected light crossing time of the whole box "
                        "(default 1/c for a 1 cm box)")
    p.add_argument("--arrival-tol", type=float, default=0.05)
    p.add_argument("--depth-tol", type=float, default=1.0e-3)
    p.add_argument("--arrival-dumps", nargs="*", default=None,
                   help="a separate, finely sampled series used for the arrival "
                        "time only (the depth series can be coarse)")
    p.add_argument("--control", default=None,
                   help="final dump of the closure = eddington control")
    p.add_argument("--json", default=None)
    p.add_argument("--map-nx", type=int, default=140)
    p.add_argument("--map-ny", type=int, default=40)
    p.add_argument("--selftest-fail", action="store_true")
    args = p.parse_args()

    args.e_in = args.arad * args.tr ** 4
    if args.t_cross is None:
        args.t_cross = 1.0 / args.c

    if args.selftest:
        dumps = selftest(args)
    else:
        if not args.dumps:
            p.error("dumps are required unless --selftest is given")
        dumps = common.load_series(args.dumps, args.bin_convert_dir,
                                   args.all_ranks)
    dumps = sorted(dumps, key=lambda d: d.time)

    if args.arrival_dumps:
        adumps = sorted(common.load_series(args.arrival_dumps,
                                           args.bin_convert_dir,
                                           args.all_ranks),
                        key=lambda d: d.time)
    else:
        adumps = dumps
    atimes = np.array([d.time for d in adumps])
    times = np.array([d.time for d in dumps])
    lit = np.array([probe(d, args.x_arrival, args.y_lit) for d in adumps])
    dark = np.array([probe(d, args.x_shadow, args.y_shadow) for d in dumps])
    litx = np.array([probe(d, args.x_shadow, args.y_lit) for d in dumps])

    thr = args.arrival_frac * args.e_in
    t_arr = float("nan")
    above = lit >= thr
    if above.any() and not above[0]:
        k = int(np.argmax(above))
        t0, t1 = atimes[k - 1], atimes[k]
        e0, e1 = lit[k - 1], lit[k]
        t_arr = t0 + (thr - e0) * (t1 - t0) / (e1 - e0)
    elif above.all():
        t_arr = atimes[0]
    t_exp = args.x_arrival * args.t_cross

    # The front is smeared over tens of cells by PLM + HLL, so a single
    # threshold crossing at one probe lags by roughly half that width.  The
    # SPEED is the invariant number: the position of the same threshold in
    # the lit lane, fitted linearly against time over the dumps in which the
    # front is inside the domain.
    xf, tf = [], []
    for d, t in zip(adumps, atimes):
        i = int(np.argmin(np.abs(np.asarray(d.x2) - args.y_lit)))
        prof = np.asarray(d.var("m1_e"))[0, i, :]
        over = prof >= thr
        if not over.any():
            continue
        k = int(np.max(np.where(over)[0]))
        if k == 0 or k >= prof.size - 1:
            continue
        # linear interpolation of the crossing between cells k and k+1
        e0, e1 = prof[k], prof[k + 1]
        w = (e0 - thr) / (e0 - e1) if e0 != e1 else 0.0
        xf.append(float(d.x1[k] + w * (d.x1[k + 1] - d.x1[k])))
        tf.append(float(t))
    speed = float("nan")
    if len(xf) >= 3:
        speed, _ = common.linfit(np.array(tf), np.array(xf))
    err_t = abs(t_arr / t_exp - 1.0) if t_exp > 0 else float("nan")

    depth = dark / np.maximum(litx, 1e-300)
    d_final = float(depth[-1])

    common.report(args, "front at x=%.3g: arrival %.6e s vs %.6e s "
                        "(%+.2f %%)" % (args.x_arrival, t_arr, t_exp,
                                        100.0 * (t_arr / t_exp - 1.0)))
    common.report(args, "front SPEED from the %.2f threshold in the lit lane "
                        "(%d dumps): %.6e cm/s = %.4f c" %
                  (args.arrival_frac, len(xf), speed, speed / args.c))
    common.report(args, "shadow depth E(%.2g,%.2g)/E(%.2g,%.2g) over time: %s"
                  % (args.x_shadow, args.y_shadow, args.x_shadow,
                     args.y_lit,
                     " ".join("%.3e" % v for v in depth)))

    d_ctrl = float("nan")
    ctrl = None
    if args.control:
        ctrl = common.load_dump(args.control, args.bin_convert_dir,
                                args.all_ranks)
        d_ctrl = (probe(ctrl, args.x_shadow, args.y_shadow)
                  / max(probe(ctrl, args.x_shadow, args.y_lit), 1e-300))
        common.report(args, "control (eddington) shadow depth = %.3e "
                            "(%.1fx the M1 value)"
                      % (d_ctrl, d_ctrl / max(d_final, 1e-300)))

    if args.json:
        import json

        def thin(d):
            arr = np.asarray(d.data["m1_e"][0])
            iy = np.linspace(0, arr.shape[0] - 1, min(args.map_ny,
                                                      arr.shape[0]))
            ix = np.linspace(0, arr.shape[1] - 1, min(args.map_nx,
                                                      arr.shape[1]))
            iy = np.unique(iy.astype(int))
            ix = np.unique(ix.astype(int))
            return ([float(v) for v in np.asarray(d.x1)[ix]],
                    [float(v) for v in np.asarray(d.x2)[iy]],
                    [[float(v) for v in row] for row in arr[np.ix_(iy, ix)]])

        mx, my, mm = thin(dumps[-1])
        yc, ec = ycut(dumps[-1], args.x_shadow)
        out = {"x": mx, "y": my, "map": mm,
               "series": {"ycut_m1": [float(v) for v in ec],
                          "ycut_y": [float(v) for v in yc],
                          "depth_t": [float(v) for v in times],
                          "depth": [float(v) for v in depth]},
               "meta": {"t_arrival": t_arr, "t_expected": t_exp,
                        "front_speed_over_c": speed / args.c,
                        "depth_final": d_final, "depth_eddington": d_ctrl,
                        "e_in": args.e_in,
                        "e_ambient": args.arad * args.t0 ** 4,
                        "time": float(times[-1])}}
        if ctrl is not None:
            cx, cy, cm = thin(ctrl)
            out["map_eddington"] = cm
            out["x_eddington"] = cx
            out["y_eddington"] = cy
            ycc, ecc = ycut(ctrl, args.x_shadow)
            out["series"]["ycut_eddington"] = [float(v) for v in ecc]
        with open(args.json, "w") as fp:
            json.dump(out, fp)

    ok = (err_t < args.arrival_tol) and (d_final < args.depth_tol)
    common.verdict(ok, "T2 shadow: arrival %.4e s vs %.4e s (%.3e, <= %.3g) "
                       "[front speed %.4f c]  depth=%.3e (<= %.3g)  "
                       "eddington control=%.3e"
                   % (t_arr, t_exp, err_t, args.arrival_tol,
                      speed / args.c, d_final, args.depth_tol, d_ctrl))


if __name__ == "__main__":
    main()
