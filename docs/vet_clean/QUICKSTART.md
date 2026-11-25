# Quick start: the clean VET radiation scheme in AthenaK

From an empty directory to a checked binary, the test suite, plots and a 3-D star.
Run the blocks in order. Each block has one sentence on what it does.

## 1. What you need

| for | what | version |
|---|---|---|
| CPU build | C++17 compiler, MPI (`mpicxx`), CMake | GCC 11 or newer (13 used), Open MPI 4.1 or 5.0, CMake 3.20 or newer |
| NVIDIA GPU | CUDA toolkit + host GCC | CUDA 12.6 or 12.9 with GCC 13 (not CUDA 12.2, not CUDA 13; see 3b) |
| AMD GPU | ROCm (`hipcc`) | ROCm 7.2 |
| tests | python3 + numpy | any recent |
| plots | matplotlib, jupyter | any recent |

On a cluster, load these as modules first (examples in section 3). The Python packages:

```bash
python3 -m pip install --user numpy matplotlib jupyter
```

## 2. Download

Clone the repository, switch to the release branch, and fetch Kokkos (a git submodule):

```bash
git clone <repo-url> athenak
cd athenak
git checkout vet-clean-release
git submodule update --init kokkos
```

All commands below run from this `athenak/` directory.

## 3. Build

Every block below makes two binaries from the same flags:

- `build_tests/src/athena`, with the built-in problem generators: for the test suite
  (`tests/vet_clean`) and the notebooks;
- `build_star/src/athena`, with `-D PROBLEM=star_vet`: for the 3-D star decks. It
  contains no built-in test problems.

`-D CMAKE_CXX_FLAGS=-ffp-contract=off` makes CPU and GPU agree to round-off; keep it.
Pick the block for your machine. On a cluster, run the GPU builds inside a compute job, not
on a login node (many sites kill long login-node processes).

### 3a. CPU (any x86 or ARM cluster or workstation)

Configure and compile with MPI and OpenMP (about 2 min per binary with 16 make jobs):

```bash
OPTS="-D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx
  -D Athena_ENABLE_MPI=ON -D Athena_ENABLE_OPENMP=ON
  -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_VET_VGD_NMAX=4096"
cmake -S . -B build_tests $OPTS
cmake -S . -B build_star $OPTS -D PROBLEM=star_vet
make -C build_tests -j 16 athena
make -C build_star -j 16 athena
```

Example modules: `module load gcc/13 openmpi/4.1 cmake/3.30`.

### 3b. NVIDIA A100, H100, H200 (CUDA)

Compile through Kokkos' `nvcc_wrapper` with MPI as the host compiler (about 6 min per
binary with 32 make jobs; set `ARCH=AMPERE80` for A100, `ARCH=HOPPER90` for H100 and H200):

```bash
ARCH=AMPERE80
export NVCC_WRAPPER_DEFAULT_COMPILER=mpicxx
OPTS="-D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_$ARCH=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_VET_VGD_NMAX=512"
cmake -S . -B build_tests $OPTS
cmake -S . -B build_star $OPTS -D PROBLEM=star_vet
make -C build_tests -j 32 athena
make -C build_star -j 32 athena
```

Example modules: `module load gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake/3.30` (A100);
on H200 use CUDA 12.9 with GCC 13. `mpicxx --showme` must name `g++`. Check that
`grep -r -- --fmad=false build_star/src/CMakeFiles/athena.dir/flags.make` finds the flag.

**Do not build with CUDA 12.2 + GCC 11 on H200.** That binary runs without errors but its
ray closure is wrong (the Eddington tensor stays at I/3 in the thin atmosphere). Use CUDA
12.9 + GCC 13. Check any new GPU binary with section 4.

### 3c. NVIDIA GH200 (aarch64, CUDA 12.9)

The same as 3b with `ARCH=HOPPER90`; on a Cray system the host compiler is `CC`:<sup>1</sup>

```bash
ARCH=HOPPER90
export NVCC_WRAPPER_DEFAULT_COMPILER=CC
OPTS="-D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_$ARCH=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_VET_VGD_NMAX=512"
cmake -S . -B build_tests $OPTS
cmake -S . -B build_star $OPTS -D PROBLEM=star_vet
make -C build_tests -j 32 athena
make -C build_star -j 32 athena
```

Example modules: `cudatoolkit` 12.9 in the Cray programming environment.

### 3d. AMD MI300A, MI250 (HIP)

Compile with `hipcc` (about 7 min per binary with 32 make jobs; set `ARCH=AMD_GFX942_APU`
for MI300A, `AMD_GFX942` for MI300X, `AMD_GFX90A` for MI250):

```bash
ARCH=AMD_GFX942_APU
OPTS="-D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=$(which hipcc)
  -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_$ARCH=On -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_VET_VGD_NMAX=512"
cmake -S . -B build_tests $OPTS
cmake -S . -B build_star $OPTS -D PROBLEM=star_vet
make -C build_tests -j 32 athena
make -C build_star -j 32 athena
```

Example modules: `module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0`.
On MI300A set `export HSA_XNACK=1` before every run, and run at most 2 MPI ranks per
MI300A.<sup>2</sup>

More options, toolchains and build times: [BUILD.md](BUILD.md).

## 4. Check the build

Run one short test, energy, mass and angular-momentum conservation in a 3-D rotating shell
(10 s on 4 cores):

```bash
tests/vet_clean/run_tests.sh --exe build_tests/src/athena --only conservation_sp
```

It must print a `PASS` line with the drifts `dE`, `dM`, `dLz` at round-off (below 1e-13):

```
PASS conservation_sp       6 s  dE 1e-15 dM 6.74e-16 dLz 1.47e-16 dLz_gr 5.92e-07
1 passed, 0 failed, 14 s
```

Then run the quick suite, 21 tests (50 s on 48 cores). By default the suite runs one test
per 4 cores of the machine; on a shared login node cap it with `--jobs` (12 jobs x 4
threads = 48 cores):

```bash
tests/vet_clean/run_tests.sh --exe build_tests/src/athena --quick --jobs 12
```

It ends with `20 passed, 0 failed`. Each line reads `PASS|FAIL  test  wall-time  measured
values`; a `FAIL` line names the failed criterion after `<-` and keeps the run directory in
`vet_clean_runs/<test>/` (its `run.log` has the details). The full suite (50 tests, 5 min on
48 cores) runs without `--quick`. What each test checks: [tests/vet_clean/README.md](../../tests/vet_clean/README.md).

On a GPU, run the same suite with the GPU `build_tests` binary, one test at a time (2 min
on one MI300A):

```bash
tests/vet_clean/run_tests.sh --exe build_tests/src/athena --threads 1 --jobs 1 --quick
```

## 5. Run a test problem by hand

A closed rotating radiating shell on a 3-D spherical-polar grid, 40 cycles in about 15 s:

```bash
mkdir -p run/shell && cd run/shell
OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=4 \
  ../../build_tests/src/athena -i ../../tests/vet_clean/decks/shell_sp.athinput
cd ../..
```

Started without `mpirun`, Open MPI binds the process to one core, so the 4 threads share it
(minutes instead of seconds); `OMPI_MCA_hwloc_base_binding_policy=none` turns that off. With
`mpirun -np 1 --bind-to none` it is not needed.

It writes history files (`sp.*.hst`) and binary dumps (`bin/*.bin`); any key can be changed
on the command line as `block/key=value`, for example `time/nlim=100`. The output files:
[OUTPUTS.md](OUTPUTS.md).

## 6. Plots in Jupyter

One notebook per test problem runs the deck and plots the result against its reference.
Name the binary in `ATHENA_EXE` (the notebooks' default is `build_vettest/src/athena`, the
binary of `run_tests.sh --build`), then run a notebook without a browser (about 1 min):

```bash
export ATHENA_EXE=$PWD/build_tests/src/athena
cd tests/vet_clean/notebooks
jupyter nbconvert --to notebook --execute --inplace shell.ipynb
cd ../../..
```

or open them interactively with `jupyter notebook` (or `jupyter lab`) in that directory.

All ten notebooks run on Python 3.6 (numpy 1.19, matplotlib 3.3); `jupyter lab` needs
JupyterLab, which `pip install jupyter` does not install on Python 3.6. How to start them
and what each one shows:
[tests/vet_clean/notebooks/README.md](../../tests/vet_clean/notebooks/README.md).

## 7. 3-D stars in spherical-polar and cubed-sphere coordinates

The star decks run with `build_star/src/athena` (section 3). The decks, the tables and the launch lines:
[STARS.md](STARS.md).

## 8. Further reading

| document | for |
|---|---|
| [RUN.md](RUN.md) | runs, production decks, restarts, step retries |
| [KEYS.md](KEYS.md) | the input keys a user sets |
| [OUTPUTS.md](OUTPUTS.md) | history columns, energy checks, `CLEAN WARNING` messages |
| [SCHEME.md](SCHEME.md) | the method and its properties |
| [BUILD.md](BUILD.md) | all build options and known-good toolchains |

---

<sup>1</sup> Block 3c is taken from the build scripts used for production on GH200; it was
not re-run for this release. Sections 2, 3a, 4, 5 and 6 and the CPU start-up check of
STARS.md were run as written from a fresh clone on a CPU cluster (GCC 13, Open MPI 4.1). Block 3b was built as written on A100 (CUDA 12.6,
GCC 13). Block 3d was built on MI300A (ROCm 7.2), and the GPU quick suite of section 4
(its 8 tests before `wedge_sym` and `floor_pay` were added) passed there on one MI300A.

<sup>2</sup> With ROCm 7.2, more than 2 ranks per MI300A on a cubed-sphere mesh hang unless
`GPU_MAX_HW_QUEUES=2` is set; the code stops with a FATAL message that says so.
