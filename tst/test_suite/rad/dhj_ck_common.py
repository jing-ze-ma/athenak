"""
Shared setup for the deep_hot_jupiter_rt correlated-k tests.

Not a test module.  It exists because those tests are unlike the rest of the suite: the
scheme lives in a USER problem generator, so they need their own binary built with
-D PROBLEM=deep_hot_jupiter_rt rather than the shared build run_test_suite.py makes, and
they need the Exo-FMS correlated-k tables, which are not in git (no upstream licence;
data/exo_fms_ck/PROVENANCE.md records where to fetch them).  Without the tables the tests
SKIP rather than fail, which is what happens in CI.
"""

# Modules
import os
import re
import subprocess
import numpy as np
import pytest

SIGMA_SB = 5.6704e-5                       # Stefan-Boltzmann, cgs

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATA = os.path.join(REPO, "data", "exo_fms_ck")
KTABLE = os.path.join(DATA, "ck", "Premixed_1x_g8_11.txt")
INPUT = os.path.join(REPO, "inputs", "mhd", "deep_hot_jupiter_rt_eos.athinput")

HAVE_TABLES = os.path.exists(KTABLE)
NO_TABLES = "Exo-FMS correlated-k tables absent; see data/exo_fms_ck/PROVENANCE.md"

# nx1 must equal the meshblock nx1 -- the two-stream RT sweeps whole radial columns -- and
# the polar boundary needs an even number of meshblocks in x3.
CK = [
    "problem/rt_ck=true",
    "problem/ck_table=" + KTABLE,
    "problem/ck_data_dir=" + DATA,
]


def mesh(nx2, nx3, mbx2, mbx3, nlim):
    """Grid and cycle-count overrides."""
    return [
        "mesh/nx1=64", "mesh/nx2=" + repr(nx2), "mesh/nx3=" + repr(nx3),
        "meshblock/nx1=64", "meshblock/nx2=" + repr(mbx2),
        "meshblock/nx3=" + repr(mbx3),
        "time/nlim=" + repr(nlim),
    ]


def build(where, flags=[]):
    """Configure and compile a deep_hot_jupiter_rt binary into `where`."""
    subprocess.run(
        ["cmake", "-S", REPO, "-B", where, "-D", "PROBLEM=deep_hot_jupiter_rt",
         "-D", "CMAKE_BUILD_TYPE=Release"] + flags,
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        ["make", "-C", where, "-j", str(os.cpu_count())],
        check=True, capture_output=True, text=True,
    )
    return os.path.join(where, "src", "athena")


def run(binary, rundir, args, ranks=0):
    """Run the problem in its own directory.  ranks = 0 means no mpirun at all."""
    os.makedirs(rundir, exist_ok=True)
    launch = ["mpirun", "-np", repr(ranks)] if ranks else []
    proc = subprocess.run(
        launch + [binary, "-i", INPUT] + args,
        cwd=rundir, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        pytest.fail(f"run failed ({ranks or 1} rank(s)):\n{proc.stdout[-4000:]}")
    return proc.stdout


def read_column(path):
    """Return (header dict, columns) from a correlated-k column dump."""
    with open(path) as f:
        head = "".join([ln for ln in f if ln.startswith("#")])
    hdr = {}
    for k in ["mu0", "icut", "T_int", "T_irr"]:
        m = re.search(k + r" = *([-0-9.eE+]+)", head)
        assert m is not None, f"'{k}' missing from the column dump header"
        hdr[k] = float(m.group(1))
    hdr["icut"] = int(hdr["icut"])
    return hdr, np.loadtxt(path)


def grab(stdout, pattern, name):
    """Pull one float out of the run's diagnostic output."""
    m = re.search(pattern, stdout)
    if m is None:
        pytest.fail(f"{name} not reported by the run; the diagnostic may have moved")
    return float(m.group(1))
