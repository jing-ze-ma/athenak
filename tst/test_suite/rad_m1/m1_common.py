"""
Shared setup for the implicit M1 regression tests (tests_m1/gates/README.md).

Not a test module.  The implicit M1 transport is exercised by the He slab of the
box_convection USER problem generator, so these tests build their own binary with
-D PROBLEM=box_convection (serial for the _cpu test, MPI for the _mpicpu ones) instead
of using the shared build of run_test_suite.py.  The binary is built once per session
and build type, and removed at the end of the session (conftest.py).

The slab needs four data files that are not in git (the He Rosseland and Planck tables
and two column profiles of bench/m1_stage2).  They are looked for in $ATHENAK_M1_DATA
(a directory holding them by name), then at their places in the checkout layout
(data/stellar_opac, ../bench); without them the tests SKIP.
"""

# Modules
import os
import shutil
import subprocess
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "he_slab_m1.athinput")
BENCH = os.path.join(REPO, "..", "bench")

# <problem> key: (file name, places to look besides $ATHENAK_M1_DATA)
FILES = {
    "ic_profile": ("ic_m1_V3edd_pgen.txt", [os.path.join(BENCH, "m1_stage2", "ic")]),
    "wb_arad_file": ("arad_V3edd.txt", [os.path.join(BENCH, "m1_stage2", "ic")]),
    "opac_table": ("rosseland_he_x0.0_z0.02.txt",
                   [os.path.join(BENCH, "hestar_fecz"),
                    os.path.join(REPO, "data", "stellar_opac")]),
    "planck_table": ("planck_he_x0.0_z0.02_ferg+tops.txt",
                     [os.path.join(REPO, "data", "stellar_opac")]),
    "m1_ic_file": ("m1_rad_ic_V3edd.txt",
                   [os.path.join(REPO, "tests_m1", "runs_3d_edd")]),
}


def _find(name, places):
    env = os.environ.get("ATHENAK_M1_DATA")
    for d in ([env] if env else []) + places:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return os.path.abspath(p)
    return None


DATA = {k: _find(n, p) for k, (n, p) in FILES.items()}
HAVE_DATA = all(DATA.values())
NO_DATA = ("He slab data files absent (set ATHENAK_M1_DATA): "
           + ", ".join(FILES[k][0] for k, v in DATA.items() if v is None))
HAVE_MPI = shutil.which("mpirun") is not None

BUILD = {False: os.path.join(REPO, "tst", "build_m1"),
         True: os.path.join(REPO, "tst", "build_m1_mpi")}

# the be levers (the default since m1-accmerge), named, and the tight solver
LEVERS = [
    "rad_m1/implicit_vimp_fold=true", "rad_m1/implicit_fast_kernels=true",
    "rad_m1/implicit_one_pass=4", "rad_m1/implicit_predictor_order=2",
]
# implicit_lin_tol 1e-12 stalls at round-off; 1e-11 converges (tests_m1/gates/gates.py)
TIGHT = ["rad_m1/implicit_tol=1.0e-12", "rad_m1/implicit_lin_tol=1.0e-11"]
VET = ["rad_m1/closure=vet_sc", "rad_m1/vet_tensor=full"]


def build(mpi):
    """Configure and compile a box_convection binary (once per session)."""
    where = BUILD[mpi]
    exe = os.path.join(where, "src", "athena")
    if os.path.exists(exe):
        return exe
    flags = ["-D", "Athena_ENABLE_MPI=ON"] if mpi else []
    subprocess.run(
        ["cmake", "-S", REPO, "-B", where, "-D", "PROBLEM=box_convection",
         "-D", "CMAKE_BUILD_TYPE=Release"] + flags,
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        ["make", "-C", where, "-j", str(os.cpu_count())],
        check=True, capture_output=True, text=True,
    )
    return exe


def run(exe, rundir, args, ranks=0, restart=None):
    """Run the slab in its own directory (emptied first unless this is a restart: the
    history files are appended to); ranks = 0 means no mpirun at all."""
    if not restart:
        shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir, exist_ok=True)
    launch = ["mpirun", "-np", repr(ranks), "--oversubscribe"] if ranks else []
    if restart:
        cmd = launch + [exe, "-r", restart, "-d", rundir] + args
    else:
        data = ["problem/%s=%s" % (k, v) for k, v in DATA.items()]
        cmd = launch + [exe, "-i", INPUT, "-d", rundir, "problem/vpert=1.0e-2"] \
            + data + args
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True, env=env)
    if proc.returncode != 0 or "Terminating on cycle limit" not in proc.stdout:
        pytest.fail(f"run failed ({ranks or 1} rank(s)):\n{proc.stdout[-4000:]}")
    return proc.stdout


def history(rundir, name="m1slab.hydro.hst"):
    """(column names, rows) of a history file."""
    names, rows = [], []
    with open(os.path.join(rundir, name)) as f:
        for ln in f:
            if ln.startswith("#"):
                if "[1]=" in ln:
                    names = [c.split("=")[1] for c in ln[1:].split()]
                continue
            rows.append([float(x) for x in ln.split()])
    return names, rows


def payload(path):
    """The bytes of a restart file after the parameter dump."""
    with open(path, "rb") as f:
        return f.read().split(b"<par_end>", 1)[-1]
