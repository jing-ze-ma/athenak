# Building AthenaK with the clean M1 radiation scheme

The clean M1 scheme is part of the normal AthenaK build: there is no extra switch. This page
gives the configurations that have been used for production and for the test suite. All of
them are plain CMake builds; out-of-source build directories are required.

## 1. Get the source

```bash
git clone <this repository> athenak
cd athenak
git checkout m1-clean-release
git submodule update --init kokkos       # Kokkos 4.6.2, pinned by the repository
```

## 2. Options that matter for this scheme

| CMake option | Value used | Why |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` turns on Kokkos bounds checks (10-50x slower). |
| `CMAKE_CXX_FLAGS` | `-ffp-contract=off` | No fused multiply-add contraction. CPU and GPU binaries then agree to round-off, and runs are reproducible across builds. On CUDA the top-level `CMakeLists.txt` also adds `--fmad=false` for the device code. Set `-D Athena_DEVICE_FMA=ON` to keep device FMA. |
| `Athena_ENABLE_MPI` | `ON` | Multi-rank runs. The CPU test suite also runs fine with MPI on. |
| `Athena_ENABLE_OPENMP` | `ON` on CPU | OpenMP threads inside a rank (the test suite uses 4 per run). |
| `Athena_M1_VGD_NMAX` | `4096` on CPU, default (432) or `512` on GPU | The largest number of ray directions of the VET closure. The star decks use 126; the test decks use up to 512. |
| `PROBLEM` | unset (built-in problem generators) or `he_star_m1` | The test suite uses only built-in generators (`rad_m1_beam`). The star decks need `-D PROBLEM=he_star_m1`. |

A binary built for one `PROBLEM` cannot run a deck of another user problem generator. The
built-in generators are always available.

## 3. CPU (workstation or cluster)

You need a C++17 compiler (GCC 11 or newer; GCC 13 is used here), MPI and CMake 3.20 or newer.

```bash
cmake -S . -B build_cpu -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx \
  -D Athena_ENABLE_MPI=ON -D Athena_ENABLE_OPENMP=ON \
  -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_M1_VGD_NMAX=4096
make -C build_cpu -j athena                # -> build_cpu/src/athena
```

The same commands are in `tests/m1_clean/run_tests.sh --build`, which builds into
`build_m1test/`.

Tested on MPCDF Raven (Intel Xeon IceLake, 72 cores per node), 2026-11-23:
`module load gcc/13 openmpi/4.1 cmake/3.30`; the configure step takes 10 s and `make -j 72`
takes 81 s.

For the star problem generator add `-D PROBLEM=he_star_m1`.

## 4. NVIDIA GPUs (CUDA)

Compile through Kokkos' `nvcc_wrapper` and name the GPU architecture:

```bash
export NVCC_WRAPPER_DEFAULT_COMPILER=mpicxx      # host compiler = the MPI C++ wrapper
cmake -S . -B build_cuda -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper \
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_<ARCH>=On \
  -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off \
  -D Athena_M1_VGD_NMAX=512 -D PROBLEM=he_star_m1
make -C build_cuda -j athena
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
must agree to a few digits. The test suite (`tests/m1_clean/run_tests.sh --exe <gpu
binary>`) is the stronger check.

Known-good module sets:

| Machine | GPU | Modules | `<ARCH>` |
|---|---|---|---|
| MPCDF Raven | A100 | `gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake/3.30` | `AMPERE80` |
| Caltech HPC (Resnick) | H200 | `gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi` (cmake from the system or a module of version >= 3.20) | `HOPPER90` |
| NCSA DeltaAI | GH200 | `cudatoolkit` 12.9 (Cray PE, host compiler `CC`: `NVCC_WRAPPER_DEFAULT_COMPILER=CC`) | `HOPPER90` |

At run time on the Caltech cluster, set `UCX_TLS=^xpmem OMPI_MCA_io=romio321`. On Raven A100
nodes, set `UCX_IB_GPU_DIRECT_RDMA=0`.

## 5. AMD GPUs (HIP, MI300A)

```bash
cmake -S . -B build_hip -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=$(which hipcc) \
  -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
  -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF \
  -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off \
  -D Athena_M1_VGD_NMAX=512 -D PROBLEM=he_star_m1
make -C build_hip -j 32 athena
```

Tested on MPCDF viper (MI300A APU) with `module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0`.
`GFX942_APU` treats host and device memory as one pool. For a discrete MI300X use
`Kokkos_ARCH_AMD_GFX942`. On ROCm 7.2, more than 2 MPI ranks per GPU on a cubed-sphere mesh
hang in the ray-sweep halo exchange unless `GPU_MAX_HW_QUEUES=2` is set; the code stops
with a FATAL message that says so.
The default UCX settings are both bitwise and fastest on viper; no extra environment is
needed.

## 6. Build time

The longest translation units are `rad_m1_vetgd.cpp` (about 12 min with nvcc, mostly one
`ptxas` process), then `rad_m1_init.cpp` and `hydro_fluxes.cpp`. With Ninja (`-G Ninja`) and
`ccache` (`-D CMAKE_CXX_COMPILER_LAUNCHER=ccache`), a rebuild after a one-file change takes
1-2 min on A100 and well under a minute on CPU. A cold CPU build takes 1.5 min on 72 cores.

## 7. After building

Run the CPU test suite (about 1.5 min with `--quick` on 16 cores):

```bash
tests/m1_clean/run_tests.sh --exe build_cpu/src/athena --quick
```

and see [RUN.md](RUN.md) for running decks.
