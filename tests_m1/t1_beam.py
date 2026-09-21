#!/usr/bin/env python3
"""T1 (beam): transverse width of a 45-degree free-streaming beam.

Design note section 8, (T1): 128^2, 45 degrees, kappa = 0,
``thick_flux = none``.  FWHM is 24 cells with the computed eigenvalues (30
with lam = +-c; Gonzalez et al. 2007).  Pass: <= 25 cells.  We also require
positivity of E and |F| <= c E (the M1 admissibility conditions of section 4).

The profile is cut PERPENDICULAR to the beam axis at mid-domain and sampled
with bilinear interpolation, so the answer does not depend on the cell
diagonal.
"""

import numpy as np

import common


def fwhm(s, prof):
    """FWHM of a single-peaked profile by linear interpolation, or None."""
    prof = np.asarray(prof, dtype=float)
    pk = int(np.argmax(prof))
    half = 0.5 * prof[pk]
    if prof[pk] <= 0.0:
        return None

    def cross(idx_range):
        prev = pk
        for i in idx_range:
            if prof[i] <= half:
                y0, y1 = prof[prev], prof[i]
                if y1 == y0:
                    return s[i]
                return s[prev] + (half - y0) * (s[i] - s[prev]) / (y1 - y0)
            prev = i
        return None

    left = cross(range(pk - 1, -1, -1))
    right = cross(range(pk + 1, prof.size))
    if left is None or right is None:
        return None
    return float(right - left)


def analyse(dump, args):
    e = dump.var("m1_e")
    if e.shape[0] != 1:
        e = e[e.shape[0] // 2:e.shape[0] // 2 + 1]
    e2 = e[0]
    x1, x2 = dump.x1, dump.x2
    dcell = 0.5 * (dump.dx1 + dump.dx2)

    ang = np.deg2rad(args.angle)
    dvec = np.array([np.cos(ang), np.sin(ang)])
    pvec = np.array([-np.sin(ang), np.cos(ang)])
    x0 = args.x0 if args.x0 is not None else x1[0] - 0.5 * dump.dx1
    y0 = args.y0 if args.y0 is not None else x2[0] - 0.5 * dump.dx2
    cx = 0.5 * (x1[0] + x1[-1])
    cy = 0.5 * (x2[0] + x2[-1])
    along = (cx - x0) * dvec[0] + (cy - y0) * dvec[1]
    px = x0 + along * dvec[0]
    py = y0 + along * dvec[1]

    half_len = args.transverse_cells * dcell
    ns = int(2 * args.transverse_cells * args.oversample) + 1
    s = np.linspace(-half_len, half_len, ns)
    prof = common.sample_2d(x1, x2, e2, px + s * pvec[0], py + s * pvec[1])
    width = fwhm(s, prof)

    emax = float(e.max())
    emin = float(e.min())
    fmax = 0.0
    if dump.has("m1_f1"):
        f1 = dump.var("m1_f1")
        f2 = dump.var("m1_f2") if dump.has("m1_f2") else np.zeros_like(f1)
        f3 = dump.var("m1_f3") if dump.has("m1_f3") else np.zeros_like(f1)
        mag = np.sqrt(f1 ** 2 + f2 ** 2 + f3 ** 2)
        sel = dump.var("m1_e") > args.f_thresh * emax
        if sel.any():
            fmax = float((mag[sel] / (args.c * dump.var("m1_e")[sel])).max())
    return width, dcell, emin, emax, fmax


def selftest(args):
    """Synthetic 45-degree Gaussian beam; --selftest-fail widens it."""
    n = 128
    x = common.centres(0.0, 1.0, n)
    dx = x[1] - x[0]
    xx, yy = np.meshgrid(x, x, indexing="xy")
    fwhm_cells = 30.0 if args.selftest_fail else 20.0
    sig = fwhm_cells * dx / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    s = (yy - xx) / np.sqrt(2.0)
    e = np.exp(-0.5 * (s / sig) ** 2)[None, :, :]
    if args.selftest_fail:
        e = e - 1e-6 * e.max()
    nhat = 1.0 / np.sqrt(2.0)
    data = dict(m1_e=e, m1_f1=args.c * e * nhat, m1_f2=args.c * e * nhat,
                m1_f3=np.zeros_like(e))
    return common.Dump(0.0, x, x, np.array([0.5]), data, "<selftest>")


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dump", nargs="?", help="2-D bin dump with m1_e, m1_f*")
    p.add_argument("--c", type=float, default=1.0,
                   help="speed of light in code units (match the athinput)")
    p.add_argument("--angle", type=float, default=45.0,
                   help="beam direction, degrees from +x1")
    p.add_argument("--x0", type=float, default=None,
                   help="beam origin x1 (default: domain lower edge)")
    p.add_argument("--y0", type=float, default=None,
                   help="beam origin x2 (default: domain lower edge)")
    p.add_argument("--max-fwhm", type=float, default=25.0,
                   help="pass threshold on the FWHM in cells")
    p.add_argument("--neg-tol", type=float, default=1e-12,
                   help="allowed min(E) as a fraction of -max(E)")
    p.add_argument("--f-tol", type=float, default=1e-10,
                   help="allowed excess of the reduced flux above 1")
    p.add_argument("--f-thresh", type=float, default=1e-10,
                   help="only cells with E > this fraction of max(E) enter f")
    p.add_argument("--transverse-cells", type=float, default=40.0,
                   help="half-length of the transverse cut, in cells")
    p.add_argument("--oversample", type=int, default=8,
                   help="sample points per cell along the cut")
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

    width, dcell, emin, emax, fmax = analyse(dump, args)
    if width is None:
        common.verdict(False, "T1 beam: no half-maximum crossing found "
                              "(profile not single-peaked?)")
    wcells = width / dcell
    common.report(args, "t=%.6g  FWHM=%.6g (%.3f cells)  min(E)=%.3e  "
                        "max(E)=%.3e  max f=%.12f"
                  % (dump.time, width, wcells, emin, emax, fmax))
    ok_w = wcells <= args.max_fwhm
    ok_e = emin >= -args.neg_tol * emax
    ok_f = fmax <= 1.0 + args.f_tol
    common.verdict(ok_w and ok_e and ok_f,
                   "T1 beam: FWHM=%.3f cells (<= %.3g)  min(E)/max(E)=%.3e "
                   "(>= %.3g)  max f=%.12f (<= 1+%.3g)"
                   % (wcells, args.max_fwhm, emin / emax if emax else 0.0,
                      -args.neg_tol, fmax, args.f_tol))


if __name__ == "__main__":
    main()
