"""CPU gates of implicit_precond = mg_gf (runs_5t_fast5box).

  python3 gates_gf.py run  <base athena> <new athena> <rundir> [arms...]
  python3 gates_gf.py eval <rundir>

1. default bitwise: the tests_m1/gates arms box_a1, box_m2, slab_a1, slab_m2 (rbgs_fwd),
   base vs new binary (cmp.py rst_bitwise); plus mg 3 (the box production preconditioner)
   base vs new, box and slab, 1 rank;
2. the same fixed point (gates.py criteria): mg_gf (K = 2, levels 3 and 1) against
   rbgs_fwd, 1 and 2 ranks, mg_gf 1 vs 2 ranks, and on the box with vet_sc (full
   tensor) mg_gf 3 against rbgs_fwd and against mg 3.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gates"))
import cmp     # noqa: E402
import gates   # noqa: E402

GF = (" rad_m1/implicit_precond=mg_gf rad_m1/implicit_mg_levels=%d"
      " rad_m1/implicit_gf_modes=%d")
MG = " rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=3"
VS = " rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
for g in ("box", "slab"):
    for a in ("a1", "m2"):
        np_, geo, ov = gates.ARMS["%s_%s" % (g, a)]
        gates.ARMS["%s_f3%s" % (g, a)] = (np_, geo, ov + GF % (3, 2))
        gates.ARMS["%s_f1%s" % (g, a)] = (np_, geo, ov + GF % (1, 2))
    np_, geo, ov = gates.ARMS["%s_a1" % g]
    gates.ARMS["%s_mg3" % g] = (np_, geo, ov + MG)
np_, geo, ov = gates.ARMS["box_a1"]
gates.ARMS["box_va1"] = (np_, geo, ov + VS)
gates.ARMS["box_vmg3"] = (np_, geo, ov + VS + MG)
gates.ARMS["box_vf3"] = (np_, geo, ov + VS + GF % (3, 2))
gates.ARMS["slab_va1"] = gates.ARMS["slab_a1"]
gates.ARMS["slab_vmg3"] = gates.ARMS["slab_mg3"]
gates.ARMS["slab_vf3"] = gates.ARMS["slab_f3a1"]
gates.PAIRS[:] = [("gf3 vs rbgs_fwd, 1 rank", "f3a1", "a1"),
                  ("gf3 vs rbgs_fwd, 2 ranks", "f3m2", "m2"),
                  ("gf3 1 vs 2 ranks", "f3a1", "f3m2"),
                  ("gf1 vs rbgs_fwd, 1 rank", "f1a1", "a1"),
                  ("gf1 vs rbgs_fwd, 2 ranks", "f1m2", "m2"),
                  ("gf1 1 vs 2 ranks", "f1a1", "f1m2"),
                  ("gf3 vs mg3, 1 rank", "f3a1", "mg3"),
                  ("vet_sc gf3 vs rbgs_fwd", "vf3", "va1"),
                  ("vet_sc gf3 vs mg3", "vf3", "vmg3"),
                  ("vet_sc mg3 vs rbgs_fwd", "vmg3", "va1")]
DEF = ["box_a1", "box_m2", "slab_a1", "slab_m2", "box_mg3", "slab_mg3"]
NEW = (["%s_%s%s" % (g, c, a) for g in ("box", "slab") for c in ("f3", "f1")
        for a in ("a1", "m2")] + ["box_va1", "box_vmg3", "box_vf3"])

if __name__ == "__main__":
    if sys.argv[1] == "run":
        base, new, rd = sys.argv[2:5]
        arms = sys.argv[5:]
        # the keys must be named in the input: copies with them in <rad_m1>
        idir = os.path.join(rd, "inp")
        os.makedirs(idir, exist_ok=True)
        for f in ("box3d_be_x", "slab2d_plm_vimp_be_x"):
            txt = open(os.path.join(gates.HERE, f + ".athinput")).read()
            txt = txt.replace("<rad_m1>\n", "<rad_m1>\nimplicit_mg_levels = 3\n"
                              "implicit_mg_halo = true\nimplicit_gf_modes = 2\n", 1)
            open(os.path.join(idir, f + ".athinput"), "w").write(txt)
        if not arms:
            gates.run(base, os.path.join(rd, "base"), DEF, idir=idir)
            gates.run(new, rd, DEF, idir=idir)
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
