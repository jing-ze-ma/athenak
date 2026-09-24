"""gates.py: pairs of M1 configurations that must agree (tests_m1/gates/README.md, 2).

  python3 gates.py run  <athena> <rundir> [arm ...]   # run the arms (all by default)
  python3 gates.py eval <rundir>                      # judge the pairs

<athena> is a CPU MPI binary built with -D PROBLEM=box_convection; the inputs are the
box3d_be_x / slab2d_plm_vimp_be_x copies next to this script (a history row every
cycle).  Short runs (box 6 cycles, slab 12): the flow amplifies any difference in the
radiation by orders of magnitude within tens of cycles (README.md, 2).  Every arm runs
twice: at the input tolerances ("loose": implicit_tol 1e-8, implicit_lin_tol 1e-10)
and at TIGHT (implicit_tol 1e-12, implicit_lin_tol 1e-11; lin_tol 1e-12 stalls at
round-off).  A pair PASSES when, measured by cmp.py,
  1. hst_cons (tight) <= 1e-10: mass, energies, E_rad totals, fluxes, every row;
  2. row1_cons SHRINKS with the tolerance: tight <= loose/10, or tight <= 1e-12 (the
     first cycle, before the flow amplifies anything: a difference that stays put
     when the solver is tightened is not solver noise);
  3. hst_dyn (tight) <= 3e-8: kinetic energies, momenta, V1max (a coarse net: V1max
     alone moves 1.5e-8 between 1 and 4 ranks in 6 box cycles);
  4. rst_cons (tight) <= 1e-9: dens, gas energy, E_rad, F_rad in the final state;
  5. no run has a NON-CONVERGED step.
rst_dyn (the momentum field) is reported only: 5e-8 of max |rho v| in every pair.
"""
import os
import re
import subprocess
import sys
import time

import cmp

MPI = "mpirun --oversubscribe --bind-to none -np"
HERE = os.path.dirname(os.path.abspath(__file__))
BOX = ("box3d_be_x mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 "
       "mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=6")
SLAB = ("slab2d_plm_vimp_be_x meshblock/nx2=8 time/nlim=12 "
        "rad_m1/closure=vet_sc rad_m1/vet_tensor=full")
# the accmerge levers (the be default since m1-accmerge) named explicitly
LEV = ("rad_m1/implicit_vimp_fold=true rad_m1/implicit_fast_kernels=true "
       "rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2")
HM = "rad_m1/implicit_halo_mpi=true"
OV = "rad_m1/implicit_halo_overlap=true"
NOV = "rad_m1/implicit_halo_overlap=false"
FC = "rad_m1/implicit_halo_ovl_faces=true"
# implicit_lin_tol = 1e-12 stalls at round-off (100 breakdowns and 20 line-Jacobi
# fallbacks in 30 box cycles); 1e-11 converges cleanly (1-6 breakdowns, no fallback)
TIGHT = "rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-11"

# arm: (ranks, geometry, overrides)
ARMS = {
    "box_a1": (1, BOX, LEV),
    "box_k1": (1, BOX, LEV + " rad_m1/implicit_krylov_dev=2"),
    "box_m2": (2, BOX, LEV + " " + HM + " " + NOV),
    "box_o2": (2, BOX, LEV + " " + HM + " " + OV),
    "box_f2": (2, BOX, LEV + " " + HM + " " + OV + " " + FC),
    "box_g2": (2, BOX, LEV + " rad_m1/implicit_halo_mpi=false"),
    "box_m4": (4, BOX, LEV + " " + HM + " " + NOV),
    "slab_a1": (1, SLAB, LEV),
    "slab_k1": (1, SLAB, LEV + " rad_m1/implicit_krylov_dev=2"),
    "slab_m2": (2, SLAB, LEV + " " + HM + " " + NOV),
    "slab_o2": (2, SLAB, LEV + " " + HM + " " + OV),
    "slab_f2": (2, SLAB, LEV + " " + HM + " " + OV + " " + FC),
    "slab_g2": (2, SLAB, LEV + " rad_m1/implicit_halo_mpi=false"),
    "slab_m4": (4, SLAB, LEV + " " + HM + " " + NOV),
}
PAIRS = [
    ("overlap on/off", "o2", "m2"),
    ("ovl_faces on/off", "f2", "o2"),
    ("krylov_dev on/off", "k1", "a1"),
    ("halo_mpi on/off", "m2", "g2"),
    ("1 vs 2 ranks", "a1", "m2"),
    ("1 vs 4 ranks", "a1", "m4"),
    ("2 vs 4 ranks", "m2", "m4"),
]


def run(exe, rdir, arms, idir=HERE, njobs=16):
    procs = []
    todo = [(a, t) for a in arms for t in ("loose", "tight")]
    while todo or procs:
        while todo and len(procs) < njobs:
            a, t = todo.pop(0)
            np_, geo, ov = ARMS[a]
            inp, *gov = geo.split()
            d = os.path.join(rdir, "%s_%s" % (a, t))
            os.makedirs(d, exist_ok=True)
            cmd = ("%s %d %s -i %s/%s.athinput -d %s problem/vpert=1.0e-2 %s %s %s"
                   % (MPI, np_, exe, idir, inp, d, " ".join(gov), ov,
                      TIGHT if t == "tight" else ""))
            log = open(os.path.join(d, "log.txt"), "w")
            env = dict(os.environ, OMP_NUM_THREADS="1")
            procs.append(subprocess.Popen("nice " + cmd, shell=True, cwd=d, stdout=log,
                                          stderr=subprocess.STDOUT, env=env))
        time.sleep(2)
        procs = [p for p in procs if p.poll() is None]


def nonconv(d):
    try:
        txt = open(os.path.join(d, "log.txt")).read()
    except OSError:
        return -1
    m = re.search(r"NON-CONVERGED=([0-9.e+-]+)", txt)
    if not m or "Terminating on cycle limit" not in txt:
        return -1
    return int(float(m.group(1)))


def evaluate(rdir):
    nfail = 0
    print("%-5s %-18s %-9s %-19s %-9s %-9s %-9s %s" % (
        "geom", "pair", "hst_cons", "row1 loose->tight", "hst_dyn", "rst_cons",
        "rst_dyn", "verdict"))
    for geo in ("box", "slab"):
        for name, x, y in PAIRS:
            ax, ay = "%s_%s" % (geo, x), "%s_%s" % (geo, y)
            dd, bad = {}, []
            for t in ("loose", "tight"):
                da = os.path.join(rdir, "%s_%s" % (ax, t))
                db = os.path.join(rdir, "%s_%s" % (ay, t))
                if not (os.path.isdir(da) and os.path.isdir(db)):
                    break
                for d in (da, db):
                    nc = nonconv(d)
                    if nc != 0:
                        bad.append("%s: %s" % (os.path.basename(d), "no run" if nc < 0
                                               else "%d NON-CONVERGED" % nc))
                dd[t] = cmp.diff(da, db)
            if len(dd) < 2:
                continue
            lo, ti = dd["loose"], dd["tight"]
            if not ti["hst_cons"] <= 1.0e-10:
                bad.append("hst_cons (%s)" % ti["hst_cons_at"])
            if not (ti["row1_cons"] <= 1.0e-12 or ti["row1_cons"] <= lo["row1_cons"]/10):
                bad.append("no shrink (%s)" % ti["row1_cons_at"])
            if not ti["hst_dyn"] <= 3.0e-8:
                bad.append("hst_dyn (%s)" % ti["hst_dyn_at"])
            if not ti["rst_cons"] <= 1.0e-9:
                bad.append("rst_cons (%s)" % ti["rst_cons_at"])
            nfail += 1 if bad else 0
            print("%-5s %-18s %-9.1e %-8.1e->%-8.1e  %-9.1e %-9.1e %-9.1e %s %s" % (
                geo, name, ti["hst_cons"], lo["row1_cons"], ti["row1_cons"],
                ti["hst_dyn"], ti["rst_cons"], ti["rst_dyn"],
                "FAIL" if bad else "PASS", "; ".join(bad)))
    print("GATES: %s" % ("PASS" if nfail == 0 else "%d FAIL" % nfail))
    return nfail


if __name__ == "__main__":
    if sys.argv[1] == "run":
        arms = sys.argv[4:] or list(ARMS)
        run(os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]), arms)
    else:
        sys.exit(1 if evaluate(sys.argv[2]) else 0)
