"""Reference pairs (line, rbgs, mg vs rbgs_fwd) with the same binary as gates_gc.py.

  python3 refgates.py run | eval   (paths below)
"""
import sys

sys.path.insert(0, "/viper/ptmp2/jinma/wt_coarse2/tests_m1/runs_5m_precond")
import gates_mg as g   # noqa: E402

gates = g.gates
rd = "/viper/ptmp2/jinma/coarse2_0925/gates"
exe = "/viper/ptmp2/jinma/coarse2_0925/bin/athena_c4_box_convection_cpu"
if sys.argv[1] == "run":
    gates.run(exe, rd, ["box_ln1", "box_rb1", "slab_ln1", "slab_rb1", "box_ga1",
                        "slab_ga1"], idir=rd + "/inp")
else:
    gates.PAIRS[:] = [p for p in gates.PAIRS if p[1] in ("ln1", "rb1", "ga1")]
    gates.evaluate(rd)
