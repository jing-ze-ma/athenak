"""CPU gates of implicit_mg_fuse (runs_5n_mgfuse/README.md).

  python3 gates_fuse.py run  <base athena> <new athena> <rundir>
  python3 gates_fuse.py eval <rundir>

1. default bitwise: gates.py arms box_a1, box_m2, slab_a1, slab_m2 (rbgs_fwd), base vs
   new binary;
2. mg fused vs unfused (levels 4, 1 and 2 ranks, loose and tight): bitwise expected;
3. mg fused vs rbgs_fwd with the gates.py criteria (runs_5m_precond/gates_mg.py pairs).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "runs_5m_precond"))
import gates_mg as g   # noqa: E402

gates, cmp = g.gates, g.cmp
for geo in ("box", "slab"):
    for a in ("a1", "m2"):
        np_, geo_, ov = gates.ARMS["%s_g%s" % (geo, a)]
        gates.ARMS["%s_u%s" % (geo, a)] = (np_, geo_, ov + " rad_m1/implicit_mg_fuse=false")
FUSE = ["box_ga1", "box_gm2", "slab_ga1", "slab_gm2"]
UNF = ["box_ua1", "box_um2", "slab_ua1", "slab_um2"]

if __name__ == "__main__":
    if sys.argv[1] == "run":
        base, new, rd = sys.argv[2:5]
        gates.run(base, os.path.join(rd, "base"), g.DEF)
        gates.run(new, rd, g.DEF)
        idir = os.path.join(rd, "inp")
        os.makedirs(idir, exist_ok=True)
        for f in ("box3d_be_x", "slab2d_plm_vimp_be_x"):
            txt = open(os.path.join(gates.HERE, f + ".athinput")).read()
            txt = txt.replace("<rad_m1>\n", "<rad_m1>\nimplicit_mg_levels = 2\n"
                              "implicit_mg_halo = true\nimplicit_mg_fuse = true\n", 1)
            open(os.path.join(idir, f + ".athinput"), "w").write(txt)
        gates.run(new, rd, FUSE + UNF, idir=idir)
    else:
        rd = sys.argv[2]
        for a in g.DEF:
            for t in ("loose", "tight"):
                d = cmp.diff(os.path.join(rd, "base", "%s_%s" % (a, t)),
                             os.path.join(rd, "%s_%s" % (a, t)))
                print("default bitwise %-8s %-5s rst_bitwise=%s hst_cons=%.1e"
                      % (a, t, d["rst_bitwise"], d["hst_cons"]))
        for f, u in zip(FUSE, UNF):
            for t in ("loose", "tight"):
                da, db = (os.path.join(rd, "%s_%s" % (x, t)) for x in (f, u))
                d = cmp.diff(da, db)
                print("fused vs unfused %-9s %-5s rst_bitwise=%s hst_cons=%.1e "
                      "nonconv %d/%d" % (f, t, d["rst_bitwise"], d["hst_cons"],
                                         gates.nonconv(da), gates.nonconv(db)))
        gates.PAIRS[:] = [p for p in gates.PAIRS if p[1] in ("ga1", "gm2")
                          and "halo off" not in p[0]]
        sys.exit(gates.evaluate(rd))
