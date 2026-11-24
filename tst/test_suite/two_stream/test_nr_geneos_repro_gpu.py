"""
Run-to-run reproducibility of the general-EOS x1 flux kernel on an accelerator.

Under a general EOS the x1 flux kernel reconstructs the derived thermodynamic variables
and then floors the interface pressure:

    PiecewiseLinearX1(member, m, k, j, il-1, iu, wder_, dl, dr);  // thread i: dl(n,i+1)
    par_for_inner(member, il, iu, [&](const int i) {
      dl(IDPR,i) = fmax(dl(IDPR,i), eos_.pfloor);                  // thread i <- dl(n,i)
      ...

The reconstruction writes dl(n,i+1) from the thread that owns i, so the floor reads a
scratch slot ANOTHER thread wrote, and successive par_for_inner loops are not implicitly
synchronised.  Without a team barrier between them, whether the floor sees the
reconstructed value or a stale one is scheduling-dependent: the run stops being
reproducible, and a few interfaces per column silently keep an unfloored pressure.  x2 and
x3 are not affected -- those reconstructions write index i from thread i.

The test therefore runs the SAME binary on the SAME input three times and requires the
histories to be bitwise identical.  It is a GPU test because the hazard needs threads of
one team running concurrently; a serial CPU build executes par_for_inner in order and
cannot expose it.

Note on sensitivity.  A race is only ever caught probabilistically, and this one is picky:
it needs the pressure floor to BIND on some interfaces but not all (where it binds on
none, or on every one, both orderings give the same answer), and it needs enough resident
teams to make the two inner loops overlap.  The configuration below -- the shipped
general-EOS deep hot Jupiter setup at its production grid -- was measured to expose it in
3 of 3 replicate pairs before the barrier was added, and in 0 of 3 after.  Smaller grids
(16x16, 32x32) and a general-EOS linear wave were both measured NOT to expose it, so do
not shrink this in the name of speed: it would still pass, but it would no longer test
anything.

The correlated-k tables are not needed -- the shipped input defaults to the grey scheme.
"""

# Modules
import os
import shutil
import pytest
import test_suite.rad.dhj_ck_common as dhj

BUILD = os.path.join(dhj.REPO, "tst", "build_geneos_repro")
CACHE = os.path.join(dhj.REPO, "tst", "build", "CMakeCache.txt")

NRUNS = 3

# Production grid.  nx1 must equal meshblock/nx1 (the two-stream RT sweeps whole radial
# columns) and the polar boundary needs an even number of meshblocks in x3.
ARGS = [
    "mesh/nx1=64", "mesh/nx2=64", "mesh/nx3=128",
    "meshblock/nx1=64", "meshblock/nx2=16", "meshblock/nx3=16",
    "time/nlim=600", "time/tlim=1.0e30",
    "output1/dt=1.0e-30",   # history every cycle, so a difference shows in the cycle
    "output2/dt=1.0e30",    # no field dumps
]

# cmake variables forwarded from the shared build, so this binary targets whatever device
# run_test_suite.py was invoked for instead of hard-coding a backend here.
FORWARD = (
    "Kokkos_ENABLE_CUDA", "Kokkos_ENABLE_HIP", "Kokkos_ENABLE_SYCL",
    "Kokkos_ENABLE_OPENMP", "CMAKE_CXX_COMPILER", "CMAKE_HIP_ARCHITECTURES",
    "CMAKE_CUDA_ARCHITECTURES", "Athena_ENABLE_MPI",
)


def device_flags():
    """Read the shared build's cache and return its device configuration as -D flags."""
    if not os.path.exists(CACHE):
        return None
    flags, on_device = [], False
    with open(CACHE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "//")) or "=" not in line:
                continue
            name, value = line.split("=", 1)
            name = name.split(":")[0]
            keep = name in FORWARD or (
                name.startswith("Kokkos_ARCH_") and value.upper() == "ON"
            )
            if not keep or value == "":
                continue
            if value.upper() in ("OFF", "FALSE", "NO"):
                continue
            if name in ("Kokkos_ENABLE_CUDA", "Kokkos_ENABLE_HIP",
                        "Kokkos_ENABLE_SYCL"):
                on_device = True
            flags += ["-D", f"{name}={value}"]
    return flags if on_device else None


@pytest.fixture(scope="module")
def binary():
    """Build a deep_hot_jupiter_rt binary for the same device as the shared build."""
    flags = device_flags()
    if flags is None:
        pytest.skip(
            "no accelerator backend in " + CACHE + "; this test needs a device build "
            "because a serial par_for_inner cannot expose the hazard"
        )
    shutil.rmtree(BUILD, ignore_errors=True)
    return dhj.build(BUILD, flags)


def history(binary, tag):
    """Run once in its own directory and return the history file's contents."""
    rundir = os.path.join(BUILD, "run_" + tag)
    shutil.rmtree(rundir, ignore_errors=True)
    dhj.run(binary, rundir, ARGS)
    hst = [f for f in os.listdir(rundir) if f.endswith(".hst")]
    assert len(hst) == 1, f"expected one history file in {rundir}, found {hst}"
    with open(os.path.join(rundir, hst[0])) as f:
        return f.read()


def test_geneos_x1_flux_is_reproducible(binary):
    """Identical runs of an identical binary must give identical answers."""
    runs = [history(binary, repr(i)) for i in range(NRUNS)]
    for i, run in enumerate(runs[1:], start=1):
        assert run == runs[0], (
            f"run {i} differs from run 0 on identical input: the general-EOS x1 flux "
            "kernel is not reproducible.  The usual cause is a missing "
            "member.team_barrier() between the derived-variable reconstruction and the "
            "interface-pressure floor in CalculateFluxes (hydro_fluxes.cpp / "
            "mhd_fluxes.cpp), where the floor reads a scratch slot another thread wrote."
        )
