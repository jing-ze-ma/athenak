"""CPU gates of implicit_precond = mg (runs_5m_precond/README.md).

  python3 gates_mg.py run  <base athena> <new athena> <rundir>
  python3 gates_mg.py eval <rundir>

1. default bitwise: the tests_m1/gates arms box_a1, box_m2, slab_a1, slab_m2 (rbgs_fwd,
   the inputs' preconditioner) with the base and the new binary: cmp.py rst_bitwise.
2. the same fixed point: mg (levels 4, halo on) against rbgs_fwd with the gates.py
   criteria (tight hst_cons <= 1e-10, row1 shrinks, hst_dyn <= 3e-8, rst_cons <= 1e-9,
   0 NON-CONVERGED), at 1 and 2 ranks, plus mg 1 vs 2 ranks.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "gates"))
import cmp     # noqa: E402
import gates   # noqa: E402

MG = " rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=4"
for g in ("box", "slab"):
    for a in ("a1", "m2"):
        np_, geo, ov = gates.ARMS["%s_%s" % (g, a)]
        gates.ARMS["%s_g%s" % (g, a)] = (np_, geo, ov + MG)
    np_, geo, ov = gates.ARMS["%s_a1" % g]
    gates.ARMS["%s_gh1" % g] = (np_, geo, ov + MG + " rad_m1/implicit_mg_halo=false")
for g in ("box", "slab"):   # reference: other valid preconditioners vs rbgs_fwd
    np_, geo, ov = gates.ARMS["%s_a1" % g]
    gates.ARMS["%s_ln1" % g] = (np_, geo, ov + " rad_m1/implicit_precond=line")
    gates.ARMS["%s_rb1" % g] = (np_, geo, ov + " rad_m1/implicit_precond=rbgs")
gates.PAIRS[:] = [("line vs rbgs_fwd (ref)", "ln1", "a1"),
                  ("rbgs vs rbgs_fwd (ref)", "rb1", "a1"),
                  ("mg vs rbgs_fwd, 1 rank", "ga1", "a1"),
                  ("mg vs rbgs_fwd, 2 ranks", "gm2", "m2"),
                  ("mg 1 vs 2 ranks", "ga1", "gm2"),
                  ("mg halo off vs rbgs_fwd", "gh1", "a1")]
DEF = ["box_a1", "box_m2", "slab_a1", "slab_m2"]
REF = ["box_ln1", "box_rb1", "slab_ln1", "slab_rb1"]
NEW = ["box_ga1", "box_gm2", "box_gh1", "slab_ga1", "slab_gm2", "slab_gh1"]

if __name__ == "__main__":
    if sys.argv[1] == "run":
        base, new, rd = sys.argv[2:5]
        gates.run(base, os.path.join(rd, "base"), DEF)
        gates.run(new, rd, DEF)
        # the mg keys must be named in the input (command-line keys the file does not
        # name are rejected): copies with them in <rad_m1>
        idir = os.path.join(rd, "inp")
        os.makedirs(idir, exist_ok=True)
        for f in ("box3d_be_x", "slab2d_plm_vimp_be_x"):
            txt = open(os.path.join(gates.HERE, f + ".athinput")).read()
            txt = txt.replace("<rad_m1>\n", "<rad_m1>\nimplicit_mg_levels = 2\n"
                              "implicit_mg_halo = true\n", 1)
            open(os.path.join(idir, f + ".athinput"), "w").write(txt)
        gates.run(new, rd, NEW, idir=idir)
    elif sys.argv[1] == "runref":
        gates.run(sys.argv[2], sys.argv[3], REF)
    else:
        rd = sys.argv[2]
        for a in DEF:
            for t in ("loose", "tight"):
                d = cmp.diff(os.path.join(rd, "base", "%s_%s" % (a, t)),
                             os.path.join(rd, "%s_%s" % (a, t)))
                print("default bitwise %-8s %-5s rst_bitwise=%s hst_cons=%.1e"
                      % (a, t, d["rst_bitwise"], d["hst_cons"]))
        sys.exit(gates.evaluate(rd))
