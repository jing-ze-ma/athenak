# Running the general-EOS hot Jupiter on GPUs

Porting notes for `deep_hot_jupiter_rt` with `eos = general`, `general_eos = table` and
`ohmic_resistivity = eos`. Companion to [`general_eos.md`](general_eos.md), which is the
parameter reference; this file is only about getting it onto an accelerator and knowing
what has and has not been checked there.

**Everything measured below was measured on CPU** (a 112-core node, 16 MPI x 7 OpenMP, the
64 x 64 x 128 spherical-polar grid) unless it says otherwise. The code now **compiles** for
`gfx942`, but **nothing in this file has been run on a GPU.** Where a statement is an
inference rather than a measurement it says so.

---

## 1. Build

Kokkos is pinned at 4.6.2 in the submodule, which knows `AMD_GFX942` (MI300A/MI300X). On
MPCDF Viper, whose `apu` partitions are MI300A, this recipe is verified to build:

```bash
git submodule update --init --recursive
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cmake -B build \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_AMD_GFX942_APU=ON \
  -DCMAKE_CXX_COMPILER=$ROCM_PATH/bin/hipcc \
  -DCMAKE_HIP_ARCHITECTURES=gfx942 \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=deep_hot_jupiter_rt
make -C build -j
```

Two things that will otherwise waste an afternoon:

- **`hipcc` is not on `PATH`** after `module load rocm/6.3`; it lives at
  `$ROCM_PATH/bin/hipcc`. Passing the bare name fails the compiler check.
- **The MPI modules are hierarchical.** `openmpi_gpu/5.0` is invisible until a compiler is
  loaded, so a bare shell configures a no-MPI build — and that build does not compile,
  because `src/bvals/physics/bfield_bcs.cpp` calls `MPI_Allreduce` unguarded by
  `MPI_PARALLEL_ENABLED`. Non-MPI builds of this branch are broken for reasons that have
  nothing to do with the GPU.

For NVIDIA, swap in `-DKokkos_ENABLE_CUDA=On -DKokkos_ARCH_<...>=On` and the
`kokkos/bin/nvcc_wrapper` compiler. Section 2 used to warn that the resistivity module
would not survive a discrete GPU; that defect is now fixed, but nothing has been *run* on
one, so treat the first NVIDIA attempt as unexplored.

`-DPROBLEM=deep_hot_jupiter_rt` is required. `deep_hot_jupiter.cpp`, the non-RT one, has no
general-EOS support and does not compile at HEAD.

---

## 2. Host pointers on the device — FIXED

This was the one real portability defect. It is now repaired; this section records what it
was, so that the pattern is recognised if it comes back.

`KOKKOS_LAMBDA` is `[=]` (`kokkos/core/src/Kokkos_Macros.hpp`), so a lambda written inside
a member function that touches a member captures **`this`** — a host pointer — and
dereferences it on the device. That was happening in four places:

| where | what it read off the host | fix |
|---|---|---|
| `current_density.hpp` | `MeshBlockPack*` in the signature of a `KOKKOS_INLINE_FUNCTION`, then `pmesh->use_spherical_polar`, `two_d`, `three_d` and six `pcoord->` arrays | takes a `CurrentDensityGeom` POD, gathered on the host by `MakeCurrentDensityGeom()` and captured by value |
| `resistivity.cpp` | `eta_b`, `use_rkg_sts`, and `max_eta` via the non-static member functions `ResistivityEOS`/`ResistivityPerna` | local copies before each `par_for`; the three `Resistivity*()` helpers are now `static` and take `max_eta` as an argument |
| `resistivity_ct.cpp`, `resistivity_update.cpp` | the RKG weights `mu` and `nu` | local copies `mu_`, `nu_` before the kernels |
| `deep_hot_jupiter_rt.cpp` | the file-scope `bool bc_outer_maxwell` | local copy `bc_outer_maxwell_` before the `usrboundaryx1_bfieldc` kernel |

The first three were invisible to the compiler: a captured `this` is legal C++ and only
misbehaves at runtime, on hardware where the device cannot resolve a host address. **On an
APU with hardware-coherent memory (MI300A) it was legal**, which is why the module worked
on Viper and nowhere else. On a discrete GPU (MI250X, any NVIDIA card) it would have
faulted or returned garbage in every resistive-EMF and resistive-timestep kernel — and it
afflicted `perna` exactly as much as `eos`.

The fourth was different: `bc_outer_maxwell` is a namespace-scope variable, not a capture,
so hipcc rejected it outright — *"reference to `__host__` variable in `__host__ __device__`
function"*. That error, not the subtle ones, is why `deep_hot_jupiter_rt` had never been
compiled for a GPU at all.

**Verification.** The change is mechanical and alters no arithmetic, so it was checked
bitwise on CPU: pre-fix and post-fix binaries produce byte-identical output for
`inputs/tests/mhd_eos_electrons.athinput`, the same with `use_rkg_sts = true`, and 20
cycles of the `deep_hot_jupiter_rt_eos` reference input (spherical-polar geometry with
super-time-stepping engaged — the path that exercises every touched kernel). The whole
binary then links for `gfx942`.

Note what this does **not** establish: it is a compile-time and a CPU-behaviour result. No
kernel here has ever executed on a GPU.

---

## 3. Meshblock decomposition

`<meshblock>/nx1` **must** equal the mesh `nx1` — the two-stream radiative transfer sweeps
whole radial columns, and the pgen aborts at startup if it does not ("only allows one
meshblock in r direction for the RT to work properly"). So the only decomposition freedom
is in theta and phi.

The tracked reference input uses `nx2 = nx3 = 8`, giving 128 blocks. That is a CPU shape.
Accelerators want a few large blocks — at least one per GPU, ideally a small multiple of
the rank count:

| mesh 64 x 64 x 128 | meshblock | blocks | verdict on a 4-GPU node |
|---|---|---|---|
| | 64 x 8 x 8 | 128 | far too many, kernel-launch bound |
| | 64 x 32 x 32 | 8 | good: 2 per GPU |
| | 64 x 64 x 64 | 2 | idles half the node |

---

## 4. Input settings that are specific to this configuration

All of these carry their reasoning in the comments of
`inputs/mhd/deep_hot_jupiter_rt_eos.athinput`; this is the summary.

- **`<units>` is required** by the general EOS. Every factor is 1 because the problem is
  already in cgs, and `mu = 1` because composition lives in the EOS.
- **`tfloor_kelvin`, not `tfloor`.** `tfloor` floors the *code* temperature, which is `p/d`
  for an ideal gas but the EOS's own `T` here — the two scales differ by `<units>/mu`, so a
  value copied across is silently wrong. See [`general_eos.md`](general_eos.md).
- **`max_eta` is also a floor on the electron fraction**, `x_e >= 230 sqrt(T)/max_eta`. At
  `1e12` that is ~1e-8, which overrides the tabulated `x_e` below ~2000 K in the upper
  atmosphere and ~3000 K in the deep gas — most of a `Teq = 2500 K` domain. `1e14` pushes
  the cap out of the way but requires `use_rkg_sts` (section 5).
- **`bc_outer_maxwell = true`** puts the outer-ghost magnetic force in as a clamped
  effective gravity. The older unbounded form scaled as `(B^2/dr)/(rho g)` — 0.03 at 3 G but
  33 at 100 G — and past O(1) it set the ghost state rather than correcting it, collapsing
  the last two active radial shells.
- **EOS table size on the device:** 281 x 451 nodes = 15 MB at `eos_dlog = 0.05` over
  `logd [-14, 0]`, `logt [1.5, 6]`. Built on the host at startup, then copied over.

---

## 5. Super-time-stepping: when it pays

`use_rkg_sts` runs the Ohmic diffusion as an RKG2 super-step. It takes
`s = ceil(-a + sqrt((a+1)^2 + tau(3+2a)))` sub-stages with `a = 0.5` and
`tau = dt_hydro/dt_diff`, **floored at `s = 3`**.

It is **more expensive per step, always.** The floor of 3 means that when diffusion does not
bind you pay three extra sub-stages for nothing. Measured cost in wall seconds per
*simulated* second (constant eta so the cap binds everywhere):

| tau | s | STS off | STS on | gain |
|---|---|---|---|---|
| 0.03 | 3 | 0.0060 | 0.0084 | **0.71x — a 41% loss** |
| 0.29 | 3 | 0.0060 | 0.0087 | **0.69x** |
| 2.87 | 4 | 0.0174 | 0.0097 | 1.79x |
| 28.7 | 11 | 0.169 | 0.0159 | 10.7x |
| 287 | 34 | 1.674 | 0.0372 | 45x |

Break-even is `tau ~ 1.5-2`. One sub-stage costs ~0.14 of a full RK2 step here, so
`cost_on/cost_off = (1 + 0.14 s)/max(1, tau)` reproduces every row to better than 10%.

**Rule of thumb: watch `dt`. If it sits at the hydro CFL value, turn STS off.**

**Accuracy.** The super-step runs *after* the complete RK2 step (`src/driver/driver.cpp`),
i.e. Lie splitting — whereas with STS off the resistive fluxes and EMFs are tasks *inside*
the RK2 stages and are second-order coupled. A/B at fixed simulated time, `eta = 1e14`:
magnetic energies agree to `<= 1.5e-5`, the decay of `1-ME` is identical to the six digits
the history file carries, and only the kinetic energies differ, by `<= 0.8%`. Good trade for
10.7x — but re-check it if the resistivity ever shapes the field you care about rather than
just damping it.

---

## 6. Runtime

- **Create the output directory first** (`mkdir -p bin` beside the input, or whatever the
  `<output>` blocks write into). AthenaK does not create it and every rank aborts with
  `MPI_ERR_NO_SUCH_FILE`.
- **Cost of the general EOS**, ideal+perna versus general+table+eos, wall per simulated
  second: 3.4x at `bbot = 3 G`, 2.9x at 10 G. Per *cycle* it is 4.4x; the table claws some
  back by giving a larger `dt` (7.64 vs 6.33 s).
  **On a GPU this ratio may well be worse:** `T(rho,e)` is a per-cell Newton/bisection
  iteration inside `ConsToPrim`, and its iteration count varies cell to cell, so it will
  diverge across a wavefront. That is an inference, not a measurement — time it on a short
  run before sizing an allocation.
- Keep double precision. The EOS root find does adapt its tolerance under
  `Athena_SINGLE_PRECISION`, but nothing in this configuration has been tested in single.

---

## 7. What is not covered by tests

- **No GPU testing of any kind.** The regression suite runs on CPU. The code now compiles
  for `gfx942` and links, which is a real check on the section-2 fix, but a clean compile
  is not a correct run: nothing here has executed a single kernel on an accelerator.
- **No regression test for either `deep_hot_jupiter_rt` boundary change.** `tst/` builds with
  the default `PROBLEM`, so a custom pgen is not reachable from the harness. The evidence
  for those changes is a parameter sweep (independence of `dfloor` from 3 to 300 G),
  recorded in the commit messages.
- The general EOS itself *is* covered: `tst/scripts/{hydro,mhd}/*_general_eos*.py` and
  `tst/scripts/mhd/mhd_eos_electrons.py`.
