"""Restart gate of implicit_precond = mg (runs_5m_precond/README.md): the CI restart test
tst/test_suite/rad_m1/test_rad_m1_restart_mpicpu.py with implicit_precond = mg,
implicit_mg_levels = 4 on top (Eddington and vet_sc), with a given binary.

  ATHENAK_M1_DATA=... python3 rst_mg.py <athena (CPU, MPI, box_convection)> <rundir>
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tst"))
import test_suite.rad_m1.m1_common as m1                       # noqa: E402
import test_suite.rad_m1.test_rad_m1_restart_mpicpu as rt      # noqa: E402

exe, rd = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
os.makedirs(rd, exist_ok=True)
inp = os.path.join(rd, "he_slab_m1_mg.athinput")
txt = open(m1.INPUT).read().replace(
    "<rad_m1>\n", "<rad_m1>\nimplicit_mg_levels = 2\nimplicit_mg_halo = true\n", 1)
open(inp, "w").write(txt)
m1.INPUT = inp
MG = ["rad_m1/implicit_precond=mg", "rad_m1/implicit_mg_levels=4"]
ok = True
for tag, extra in (("edd", []), ("vet", m1.VET)):
    a, b = os.path.join(rd, tag + "_a"), os.path.join(rd, tag + "_b")
    args = rt.ARGS + MG + extra
    m1.run(exe, a, args + ["time/nlim=%d" % rt.N], ranks=2)
    m1.run(exe, b, args + ["time/nlim=%d" % (rt.N // 2)], ranks=2)
    rst = sorted(glob.glob(os.path.join(b, "rst", "*.rst")))[-1]
    m1.run(exe, b, ["time/nlim=%d" % rt.N], ranks=2, restart=rst)
    ra = sorted(glob.glob(os.path.join(a, "rst", "*.rst")))[-1]
    rb = sorted(glob.glob(os.path.join(b, "rst", "*.rst")))[-1]
    same = m1.payload(ra) == m1.payload(rb)
    ok = ok and same
    print("%s: %d cycles vs %d + restart + %d: final restart data %s"
          % (tag, rt.N, rt.N // 2, rt.N // 2, "BITWISE" if same else "DIFFER"))
print("RESTART GATE:", "PASS" if ok else "FAIL")
