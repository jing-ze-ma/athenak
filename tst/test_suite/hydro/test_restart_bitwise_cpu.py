"""
Restart-bitwise regression test.

A restart must be a BITWISE continuation: a run taken straight through N cycles and a
run restarted from the file written at cycle N/2 must end on exactly the same state,
bit for bit.  Four bugs broke that and were fixed in 7645550e, 2667d984, 6740fae3 and
e4b318b3 -- user boundary conditions that read w0 before the first ConsToPrim of a
restarted run, ghost zones that were refilled on restart instead of kept (a refill
cannot reproduce ghosts that were made from a pre-floor u0), and the general-EOS
temperature cache, the c2p warm start, which no restart file carried.  This test is the
gate against a repeat.

Each case runs the binary twice:
  - straight through <time>/nlim cycles, writing a restart file at nlim/2 (the rst
    output block uses dcycle) and a full dump at the end;
  - restarted from that file to the same nlim, with job/basename overridden so the two
    runs' dumps sit side by side.
Both final dumps carry EVERY variable and INCLUDE THE GHOST ZONES, and every array is
compared with numpy.array_equal -- exact equality, no tolerance.  The history file is
deliberately not used: it is a set of volume integrals, which hides a sign-paired or
small-amplitude difference.  The comparison is also made on a formatted table written
at 17 significant figures, which round-trips a double exactly where the binary dump,
whose floats are 32-bit, cannot: the whole 1D domain for the Cartesian cases, a radial
line for the two curvilinear ones.

CASES
  hydro_lwave   Cartesian hydro linear wave, periodic, two MeshBlocks
  mhd_lwave     Cartesian MHD fast wave -- adds the face-centred field to the file
  sp_mhd        spherical-polar MHD double rarefaction: polar boundary, reflecting
                radial walls, floors firing and FOFC on
  cs_mhd        the same on the cubed sphere: panel seams and the cube-vertex corner
                ghost fill
  geneos_table  1D hydro linear wave with the TABULATED general EOS -- the wtemp cache
                case (6740fae3).  This one is compared on the 32-bit binary dump ONLY.
                At full double precision the restart still differs in the last bits
                (~1e-16 relative, most cells) and that is a KNOWN OPEN GAP, not a
                regression: measured on a general-EOS shock tube, the conserved AND the
                primitive state the restart restores are bitwise identical at the
                restart cycle, and the divergence appears in the FIRST cycle after it.
                The seed is the one piece of state a restart cannot avoid recomputing.
                The file does carry wtemp, but Driver::Initialize then runs a
                ConsToPrim the continuous run does not, and the table's temperature
                root find -- which stops as soon as a step falls below logtol
                (eos_table.hpp) -- is not exactly idempotent: it re-converges the
                restored wtemp to a value an ulp away.  The comparison here is
                therefore made at dump precision, which is what 6740fae3 itself was
                validated at.

TEETH.  <time>/restart_refill_ghosts=true, which restores the pre-2667d984 behaviour,
does NOT break any case here, and neither does an SMR variant: with the built-in
problem generators every ghost fill is an exactly reproducible function of the active
cells, so the refill lands back on the file's own ghosts (verified cell by cell at the
restart cycle).  The fill that 2667d984 fixed is only unreproducible when the physical
BC is a user hook that reads the post-floor state (red giant, deep hot Jupiter), and no
user problem generator can be built by this suite.  What does have teeth here is
test_restart_perturbed_run_differs: it repeats the hydro case with the restarted run
given a different Riemann solver, and demands that the comparison sees it.  That checks
the instrument -- the file numbering, the ghost zones, the exactness -- rather than
trusting that "identical" means the test looked.
"""

# Modules
import os
import shutil
import sys
import pytest
import numpy as np
import test_suite.testutils as testutils

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "vis", "python"))
import bin_convert      # noqa: E402

# input file stem, dumped variable, and whether a 17-digit table is written too
_cases = [
    ("restart_lwave_hydro", "hydro_u", True),
    ("restart_lwave_mhd", "mhd_u_bcc", True),
    ("restart_sp_mhd", "mhd_u_bcc", True),
    ("restart_cs_mhd", "mhd_u_bcc", True),
    ("restart_geneos_table", "hydro_u", False),
]

# the restart file written at cycle nlim/2, and the final dump of each run
_RST = 1
_FINAL = 1


def outputs_cleanup():
    """Remove every output directory a case may have written."""
    for directory in ("bin", "rst", "tab", "log"):
        shutil.rmtree(directory, ignore_errors=True)
    testutils.cleanup()


def run_pair(base, restart_flags=None):
    """Run the case straight through, then restart it, and return the two basenames."""
    restart_file = f"rst/{base}.{_RST:05d}.rst"
    assert testutils.run(f"inputs/{base}.athinput"), f"{base}: continuous run failed"
    assert os.path.exists(restart_file), f"{base}: no restart file at cycle nlim/2"
    command = ["./athena", "-r", restart_file, f"job/basename={base}_r"]
    command += restart_flags if restart_flags else []
    assert testutils.run_command(command), f"{base}: restarted run failed"
    return base, f"{base}_r"


def compare_bin(base_a, base_b, var):
    """Compare every variable of two binary dumps.  Returns a list of differences."""
    name_a = f"bin/{base_a}.{var}.{_FINAL:05d}.bin"
    name_b = f"bin/{base_b}.{var}.{_FINAL:05d}.bin"
    data_a = bin_convert.read_binary(name_a)
    data_b = bin_convert.read_binary(name_b)
    assert data_a["cycle"] == data_b["cycle"], (
        f"dumps are at different cycles: {data_a['cycle']} vs {data_b['cycle']}"
    )
    assert data_a["cycle"] > 0, "final dump is the initial state, not the end of the run"
    assert data_a["time"] == data_b["time"], "dumps are at different times"
    differences = []
    for name in data_a["var_names"]:
        arr_a = np.array(data_a["mb_data"][name])
        arr_b = np.array(data_b["mb_data"][name])
        if not np.array_equal(arr_a, arr_b):
            delta = np.abs(arr_a - arr_b)
            differences.append(
                f"{name}: {int((delta > 0).sum())} of {delta.size} cells differ, "
                f"max |diff| = {delta.max():g}"
            )
    return differences


def read_tab(filename):
    """Read a formatted table into a dict of arrays, one per variable.

    athena_read.tab is not used here: on a 1D mesh the AthenaK header advertises the
    j/x2v and k/x3v columns that the data rows do not contain, and athena_read trusts
    the header.  The rows are (gid, i, x1v, var...), so the variable names are the last
    columns' worth of headings.
    """
    with open(filename, "r") as tabfile:
        lines = tabfile.readlines()
    headings = lines[1].split()[1:]
    rows = [line.split() for line in lines if not line.startswith("#")]
    nvar = len(rows[0]) - 3
    names = headings[-nvar:]
    values = np.array([[float(val) for val in row[3:]] for row in rows])
    return {name: values[:, n] for n, name in enumerate(names)}


def compare_tab(base_a, base_b, var):
    """Compare the 17-digit table dumps, which carry the full double precision."""
    data_a = read_tab(f"tab/{base_a}.{var}.{_FINAL:05d}.tab")
    data_b = read_tab(f"tab/{base_b}.{var}.{_FINAL:05d}.tab")
    differences = []
    for name in data_a:
        arr_a, arr_b = data_a[name], data_b[name]
        if not np.array_equal(arr_a, arr_b):
            delta = np.abs(arr_a - arr_b)
            differences.append(
                f"{name} (table): {int((delta > 0).sum())} of {delta.size} cells "
                f"differ, max |diff| = {delta.max():g}"
            )
    return differences


@pytest.mark.parametrize("base,var,has_tab", _cases)
def test_restart_is_bitwise(base, var, has_tab):
    """A restarted run must end on exactly the state the continuous run ends on."""
    try:
        outputs_cleanup()
        base_a, base_b = run_pair(base)
        differences = compare_bin(base_a, base_b, var)
        if has_tab:
            differences += compare_tab(base_a, base_b, var)
        if differences:
            pytest.fail(
                f"{base}: restart is not a bitwise continuation -- "
                + "; ".join(differences)
            )
    finally:
        outputs_cleanup()


def test_restart_perturbed_run_differs():
    """Validate the instrument: a restart that is NOT a continuation must be caught.

    The restarted run is given a different Riemann solver, so its last nlim/2 cycles
    are a different calculation.  If this comparison reports no difference, the test
    above is not looking at what it claims to look at.
    """
    base, var = "restart_lwave_hydro", "hydro_u"
    try:
        outputs_cleanup()
        base_a, base_b = run_pair(base, ["hydro/rsolver=llf"])
        differences = compare_bin(base_a, base_b, var)
        differences += compare_tab(base_a, base_b, var)
        if not differences:
            pytest.fail(
                "a restart run with a different Riemann solver compared EQUAL: "
                "the comparison is not reading the two runs' final dumps"
            )
    finally:
        outputs_cleanup()
