# Running the general-EOS hot Jupiter on GPUs

Porting notes for `deep_hot_jupiter_rt` with `eos = general`, `general_eos = table` and
`ohmic_resistivity = eos`. Companion to [`general_eos.md`](general_eos.md), which is the
parameter reference; this file is only about getting it onto an accelerator and knowing
what has and has not been checked there.

**Everything below was measured on CPU** (a 112-core node, 16 MPI x 7 OpenMP, the
64 x 64 x 128 spherical-polar grid) unless it says otherwise. **Nothing in this file has
been run on a GPU.** Where a statement is an inference rather than a measurement it says
so.

---

## 1. Build

Kokkos is pinned at 4.6.2 in the submodule, which knows `AMD_GFX942` (MI300A/MI300X):

```bash
git submodule update --init --recursive
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=hipcc \
  -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX942=ON \
  -DAthena_ENABLE_MPI=ON -DPROBLEM=deep_hot_jupiter_rt
make -j
```

For NVIDIA, swap in `-DKokkos_ENABLE_CUDA=On -DKokkos_ARCH_<...>=On` and the
`kokkos/bin/nvcc_wrapper` compiler — **but read section 2 first, because the resistivity
module will not work on a discrete GPU.**

`-DPROBLEM=deep_hot_jupiter_rt` is required. `deep_hot_jupiter.cpp`, the non-RT one, has no
general-EOS support and does not compile at HEAD.

---

## 2. The resistivity module dereferences host pointers on the device

This is the one real portability defect, and it is worth understanding before choosing
hardware.

`KOKKOS_LAMBDA` is `[=]` (`kokkos/core/src/Kokkos_Macros.hpp`), so a lambda written inside
a member function that touches a member captures **`this`** — a host pointer — and
dereferences it on the device. `src/diffusion/resistivity.cpp` does this in roughly a dozen
kernels, via `eta_b`, `use_rkg_sts` and `pmy_pack`.

Worse, it is not confined to the lambdas. `CurrentDensity()` in
`src/diffusion/current_density.hpp` is a `KOKKOS_INLINE_FUNCTION` that takes
`MeshBlockPack *pmy_pack` and reads `pmy_pack->pmesh->use_spherical_polar` and six
`pmy_pack->pcoord->...` arrays **on the device**. The host pointer is in the signature of
the shared helper.

- **On an APU with hardware-coherent memory (MI300A), this is legal** and the code runs.
  That is the only reason it works today.
- **On a discrete GPU (MI250X, any NVIDIA card) it is not.** Expect a fault or garbage in
  every kernel that computes a resistive EMF or the resistive timestep. This affects
  `perna` exactly as much as `eos`; it is not specific to the general EOS.

**The fix, if it is ever needed,** is contained: `current_density.hpp` is included by
`resistivity.cpp` and nothing else. Pass the three mesh flags and six coordinate Views
explicitly (a small POD built on the host and captured by value) instead of `pmy_pack`, and
take local copies of `eta_b` / `use_rkg_sts` before each `par_for`. It is mechanical and
changes no arithmetic, so it can be verified bitwise against a CPU run. Alternatively
`KOKKOS_CLASS_LAMBDA` (`[=, *this]`) fixes the lambda captures but not `CurrentDensity`.

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

- **No GPU testing of any kind.** The regression suite runs on CPU.
- **No regression test for either `deep_hot_jupiter_rt` boundary change.** `tst/` builds with
  the default `PROBLEM`, so a custom pgen is not reachable from the harness. The evidence
  for those changes is a parameter sweep (independence of `dfloor` from 3 to 300 G),
  recorded in the commit messages.
- The general EOS itself *is* covered: `tst/scripts/{hydro,mhd}/*_general_eos*.py` and
  `tst/scripts/mhd/mhd_eos_electrons.py`.
