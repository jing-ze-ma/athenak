"""CPU gates of implicit_precond = mg_gc (runs_5p_coarse2/README.md).

  python3 gates_gc.py run  <base athena> <new athena> <rundir> [arms...]
  python3 gates_gc.py eval <rundir>

1. default bitwise: the tests_m1/gates arms box_a1, box_m2, slab_a1, slab_m2 (rbgs_fwd),
   base vs new binary (cmp.py rst_bitwise);
2. the same fixed point: mg_gc (levels 1 and 3; 2 x 2 bands box, 2 bands slab) against
   rbgs_fwd with the gates.py criteria, 1 and 2 ranks, plus mg_gc 1 vs 2 ranks; and
   implicit_bcg_rho_direct (with rbgs_fwd and with mg_gc 3) against rbgs_fwd.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gates"))
import cmp     # noqa: E402
import gates   # noqa: E402

GC = " rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=%d"
RD = " rad_m1/implicit_bcg_rho_direct=true"
for g in ("box", "slab"):
    bands = (" rad_m1/implicit_gc_bands2=2 rad_m1/implicit_gc_bands3=2" if g == "box"
             else " rad_m1/implicit_gc_bands2=2")
    for a in ("a1", "m2"):
        np_, geo, ov = gates.ARMS["%s_%s" % (g, a)]
        gates.ARMS["%s_c1%s" % (g, a)] = (np_, geo, ov + GC % 1)
        gates.ARMS["%s_c3%s" % (g, a)] = (np_, geo, ov + GC % 3)
        gates.ARMS["%s_cb%s" % (g, a)] = (np_, geo, ov + GC % 1 + bands)
    np_, geo, ov = gates.ARMS["%s_a1" % g]
    gates.ARMS["%s_c3d1" % g] = (np_, geo, ov + GC % 3 + RD)
    gates.ARMS["%s_rbd1" % g] = (np_, geo, ov + RD)
gates.PAIRS[:] = [("gc1 vs rbgs_fwd, 1 rank", "c1a1", "a1"),
                  ("gc1 vs rbgs_fwd, 2 ranks", "c1m2", "m2"),
                  ("gc1 1 vs 2 ranks", "c1a1", "c1m2"),
                  ("gc3 vs rbgs_fwd, 1 rank", "c3a1", "a1"),
                  ("gc3 vs rbgs_fwd, 2 ranks", "c3m2", "m2"),
                  ("gc3 1 vs 2 ranks", "c3a1", "c3m2"),
                  ("gcb vs rbgs_fwd, 1 rank", "cba1", "a1"),
                  ("gcb vs rbgs_fwd, 2 ranks", "cbm2", "m2"),
                  ("gcb 1 vs 2 ranks", "cba1", "cbm2"),
                  ("gc3 + rho_direct vs rbgs_fwd", "c3d1", "a1"),
                  ("rbgs_fwd + rho_direct vs rbgs_fwd", "rbd1", "a1")]
DEF = ["box_a1", "box_m2", "slab_a1", "slab_m2"]
NEW = (["%s_%s%s" % (g, c, a) for g in ("box", "slab") for c in ("c1", "c3", "cb")
        for a in ("a1", "m2")] + ["%s_%s" % (g, c) for g in ("box", "slab")
                                  for c in ("c3d1", "rbd1")])

if __name__ == "__main__":
    if sys.argv[1] == "run":
        base, new, rd = sys.argv[2:5]
        arms = sys.argv[5:]
        if not arms:
            gates.run(base, os.path.join(rd, "base"), DEF)
            gates.run(new, rd, DEF)
        # the keys must be named in the input: copies with them in <rad_m1>
        idir = os.path.join(rd, "inp")
        os.makedirs(idir, exist_ok=True)
        for f in ("box3d_be_x", "slab2d_plm_vimp_be_x"):
            txt = open(os.path.join(gates.HERE, f + ".athinput")).read()
            txt = txt.replace("<rad_m1>\n", "<rad_m1>\nimplicit_mg_levels = 1\n"
                              "implicit_mg_halo = true\nimplicit_gc_bands2 = 1\n"
                              "implicit_gc_bands3 = 1\n"
                              "implicit_bcg_rho_direct = false\n", 1)
            open(os.path.join(idir, f + ".athinput"), "w").write(txt)
        gates.run(new, rd, arms or NEW, idir=idir)
    else:
        rd = sys.argv[2]
        for a in DEF:
            for t in ("loose", "tight"):
                d = cmp.diff(os.path.join(rd, "base", "%s_%s" % (a, t)),
                             os.path.join(rd, "%s_%s" % (a, t)))
                print("default bitwise %-8s %-5s rst_bitwise=%s hst_cons=%.1e"
                      % (a, t, d["rst_bitwise"], d["hst_cons"]))
        sys.exit(gates.evaluate(rd))
