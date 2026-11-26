# Building AthenaK with the clean VET radiation scheme

The clean VET scheme is part of the normal AthenaK build: there is no extra switch. This page
gives the configurations that have been used for production and for the test suite. All of
them are plain CMake builds; out-of-source build directories are required. For the shortest
path from download to a checked binary see [QUICKSTART.md](QUICKSTART.md).

## 1. Get the source

```bash
git clone <this repository> athenak
cd athenak
git checkout vet-clean-release
git submodule update --init kokkos       # Kokkos 4.6.2, pinned by the repository
```

On a cluster, build inside a compute job (`srun`/`salloc`), not on a login node: many sites
limit login-node CPU time (NCSA DeltaAI kills processes after 30 CPU-minutes), and one file
of a CUDA build (`rays/sc_sph_sweep.cpp`) takes about 12 min of nvcc alone.

## 2. Options that matter for this scheme

| CMake option | Value used | Why |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` turns on Kokkos bounds checks (10-50x slower). |
| `CMAKE_CXX_FLAGS` | `-ffp-contract=off` | No fused multiply-add contraction. CPU and GPU binaries then agree to round-off, and runs are reproducible across builds. On CUDA the top-level `CMakeLists.txt` also adds `--fmad=false` for the device code. Set `-D Athena_DEVICE_FMA=ON` to keep device FMA. |
| `Athena_ENABLE_MPI` | `ON` | Multi-rank runs. The CPU test suite also runs fine with MPI on. |
| `Athena_ENABLE_OPENMP` | `ON` on CPU | OpenMP threads inside a rank (the test suite uses 4 per run). |
| `Athena_VET_VGD_NMAX` | `4096` on CPU, default (432) or `512` on GPU | The largest number of ray directions of the VET closure. The star decks use 126; the test decks use up to 512. |
| `PROBLEM` | unset (built-in problem generators) or `star_vet` | The test suite uses only built-in generators (`rad_vet_beam`). The star decks need `-D PROBLEM=star_vet`. |

A binary built with `-D PROBLEM=<name>` runs only that problem generator: a
`star_vet` binary cannot run the test decks. Build the test binary without `-D PROBLEM`
(the GPU blocks below name `star_vet`; drop it, in a second build directory, for the
test suite on a GPU).

## 3. CPU (workstation or cluster)

You need a C++17 compiler (GCC 11 or newer; GCC 13 is used here), MPI and CMake 3.20 or newer.

```bash
cmake -S . -B build_cpu -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx \
  -D Athena_ENABLE_MPI=ON -D Athena_ENABLE_OPENMP=ON \
  -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_VET_VGD_NMAX=4096
make -C build_cpu -j 16 athena             # -> build_cpu/src/athena
```

The same commands are in `tests/vet_clean/run_tests.sh --build`, which builds into
`build_vettest/`.

Tested on an Intel Xeon IceLake cluster (72 cores per node) with
`module load gcc/13 openmpi/4.1 cmake/3.30`; the configure step takes 10 s and `make -j 72`
takes 81 s.

For the star problem generator add `-D PROBLEM=star_vet`.

## 4. NVIDIA GPUs (CUDA)

Compile through Kokkos' `nvcc_wrapper` and name the GPU architecture:

```bash
export NVCC_WRAPPER_DEFAULT_COMPILER=mpicxx      # host compiler = the MPI C++ wrapper
cmake -S . -B build_cuda -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper \
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_<ARCH>=On \
  -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off \
  -D Athena_VET_VGD_NMAX=512 -D PROBLEM=star_vet
make -C build_cuda -j 32 athena
```

`<ARCH>` is `AMPERE80` for A100, `HOPPER90` for H100 / H200 / GH200.
`Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF` turns off the stream-ordered allocator, as in every
production build (the same for HIP below). Check that `--fmad=false` appears in the build commands
(`grep fmad build_cuda/build.ninja` or the make log); it is what makes GPU and CPU agree.

**Toolchain warning (H200).** Do not use CUDA 12.2 with g++ 11 on H200. The binary builds and runs, but the
VET ray closure comes out wrong without any error message: the Eddington tensor stays at
I/3 in the optically thin atmosphere, so the photosphere never free-streams. Use CUDA 12.9 with g++ 13. CUDA 12.6
(A100) and 12.9 (A100, H200, GH200) with g++ 13 are known to be correct. CUDA 13 is not
supported by Kokkos 4.6.2. To check a new GPU binary, run one of the star decks for 2 cycles
and compare the start-up line `clean_init steady pass 1` with a CPU run of the same deck; it
must agree to a few digits. The test suite (`tests/vet_clean/run_tests.sh --exe <gpu
binary>`) is the stronger check.

Known-good module sets:

| Machine | GPU | Modules | `<ARCH>` |
|---|---|---|---|
| MPCDF Raven | A100 | `gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake/3.30` | `AMPERE80` |
| Caltech HPC (Resnick) | H200 | `gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi` (cmake from the system or a module of version >= 3.20) | `HOPPER90` |
| NCSA DeltaAI | GH200 | the default modules (`PrgEnv-gnu`, `gcc-native/14`, `cudatoolkit/25.5_12.9`, `cray-mpich`; host compiler `CC` = g++ 14: `NVCC_WRAPPER_DEFAULT_COMPILER=CC`; system cmake 3.28; configure checked for this release, build not re-run) | `HOPPER90` |

At run time on the H200 cluster above, set `UCX_TLS=^xpmem OMPI_MCA_io=romio321`; on more than one
H200 node also `UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln`. Without them the first message
between nodes fails with `ib_mlx5 ... Local protection error` (RDMA_WRITE). With `put_ppln`, two nodes of
4 H200 were bitwise from run to run. On the A100 cluster, set `UCX_IB_GPU_DIRECT_RDMA=0`.

## 5. AMD GPUs (HIP, MI300A)

```bash
cmake -S . -B build_hip -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=$(which hipcc) \
  -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
  -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF \
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off \
  -D Athena_VET_VGD_NMAX=512 -D PROBLEM=star_vet
make -C build_hip -j 32 athena
```

Tested on MI300A APUs with `module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0`.
`GFX942_APU` treats host and device memory as one pool. For a discrete MI300X use
`Kokkos_ARCH_AMD_GFX942`, for MI250 / MI250X `Kokkos_ARCH_AMD_GFX90A` (not tested). On ROCm 7.2, more than 2 MPI ranks per GPU on a cubed-sphere mesh
hang in the ray-sweep halo exchange unless `GPU_MAX_HW_QUEUES=2` is set; the code stops
with a FATAL message that says so.
With one MPI rank per APU the default UCX settings are bitwise and fastest on MI300A; no
extra environment is needed.

**Known issue: two MPI ranks sharing one APU.** When two ranks run on the same MI300A and the
APU is heavily oversubscribed (several jobs on it at once), a run can, rarely, give a wrong
result without any warning: the dumps of a 2-rank run then differ from those of a 1-rank run
from the first step on, in the last digits or by much more. Two kinds were seen in short 2-rank test runs (2 cycles,
the closed MHD shell of the test suite on spherical-polar):
- the UCX `rocm_ipc` transport between the two ranks: about 1 in 30 to 70 runs;
- a rarer kind, about 0.3 % of runs, seen only with `vet_gd_band_ring = true` (the default
  with more than one rank on a 3-D spherical-polar band); 0 in about 2400 runs with it off.

Not seen with one rank per APU (0 of 388 short runs, 0 of 94 runs of 500 cycles), on NVIDIA
H200 (0 of 3643) or on CPUs (0 of 2240). The previous release shows the first kind too. For runs
that share an APU between ranks, `UCX_TLS=^ud,tcp,rocm_ipc` removes the main kind at about
8.6 % more time per cycle. The cause of the rarer kind is not known yet.

## 6. Build time

The longest translation units are `rays/sc_sph_sweep.cpp` (about 12 min with nvcc, mostly one
`ptxas` process), then `rad_vet_params_implicit.cpp` and `hydro_fluxes.cpp`. With Ninja (`-G Ninja`) and
`ccache` (`-D CMAKE_CXX_COMPILER_LAUNCHER=ccache`), a rebuild after a one-file change takes
1-2 min on A100 and well under a minute on CPU. A cold CPU build takes 1.5 min on 72 cores.

## 7. After building

Run the CPU test suite (about 1 min with `--quick` on 48 cores):

```bash
tests/vet_clean/run_tests.sh --exe build_cpu/src/athena --quick
```

and see [RUN.md](RUN.md) for running decks.
