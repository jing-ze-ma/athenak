---
name: gpu-this-capture-device-lambda
description: FIXED af539941 + e2e1237a, and VERIFIED BY RUNNING on apudev: cubed-sphere hydro AND MHD now work on GPU and match CPU. Records the constant-memory closure limit that a 6-View POD blew (crashing MHD), and that XNACK=0 does NOT discriminate on an APU
metadata:
  type: project
---

**FIXED 2026-08-29, af539941 + e2e1237a, and VERIFIED BY ACTUALLY RUNNING** on
`apudev` (MI300A, `vipa1001`). A clean `gfx942` build gives 22 warnings -> 0 and 0
errors; cubed-sphere hydro AND MHD run on GPU and match CPU; the CPU answers are
bit-identical to 17 digits, so nothing moved.

## RESULTS ON THE GPU

| | |
|---|---|
| hydro vs CPU | agrees to **4.6e-14** on every `.hst` column; 1.2-1.4e7 zone-cycles/s |
| MHD vs CPU | convergence errors **identical to all 9 printed digits** (L1_B 5.86663501e-05) |
| XNACK 1 vs 0 | **bit-identical** for both |

## TWO MISTAKES I MADE -- both worth not repeating

**1. An INCREMENTAL build reported "0 warnings" because it never recompiled the
files that had them.** `make` in an already-built dir re-emits nothing for
untouched TUs. It hid six live sites (`GridPiecewiseLinearDerX1`,
`WbPiecewiseLinearDerX1/2/3`, `WbStaticPiecewiseLinearDerX1/2/3`).
**Always count warnings from a CLEAN build dir.**

**2. `HSA_XNACK=0` does NOT prove independence from unified memory on an APU.**
The PRE-refactor binary also passes with XNACK=0, because on MI300A host memory is
physically the same memory and is addressable from the device either way. I claimed
the XNACK=0 run validated the fix; it does not discriminate at all. The real
validation is the compile-time warning count from a clean build, plus the fact that
the code no longer names host pointers in device code.

## THE CRASH af539941 CAUSED, and the rule it teaches

`af539941` handed every flux kernel a `GnomonicTrig` with all SIX trig Views. MHD
then died on cycle 1 with
`HSA_STATUS_ERROR_MEMORY_APERTURE_VIOLATION`. **Everything a kernel captures rides
in the constant-memory buffer of `hip_parallel_launch_constant_memory`, which is
size-limited**, and `MHD::CalculateFluxes` already captures a great many Views.
Hydro survived only because its closure is smaller. `e2e1237a` shrinks the POD to
the ONE (sin, cos) pair each sweep actually needs. **Rule: adding captures to an
already-fat kernel can break it at LAUNCH, with a fault that looks like a wild
pointer. Give a kernel only what it uses.**

How it was localised, in order: `AMD_SERIALIZE_KERNEL=3 AMD_LOG_LEVEL=4` named the
dispatch (`ShaderName` demangles to `MHD::CalculateFluxes<MHD_RSolver 3>`, league
1944 = 6 blocks x 18 x 18 = the x1 sweep); Kokkos bounds checking found NOTHING,
which ruled out View indexing; a Debug build did not crash, which said
optimizer/launch rather than logic; and building **806a9e6f + only the `sync_host`
change** proved it was the refactor and not pre-existing. **That last step is the
one I nearly skipped** -- I had told the user it was pre-existing on the strength of
"it still crashes with the seam exchange disabled", which ruled out the wrong commit.

## The pattern

`CLAUDE.md` says never dereference `this` or a host-side pointer inside a device lambda.
The code does it anyway, and **hipcc says so explicitly**:

```
warning: capture host side class data member by this pointer in device or host device
lambda function may result in invalid memory access if this pointer is not accessible
on device side [-Wgpu-maybe-wrong-side]
```

`pmy_pack` is a class member, so `KOKKOS_LAMBDA` (`[=]`) captures `this` implicitly and
`pmy_pack->pmesh->use_cubed_sphere` becomes a host-pointer load in device code.

**22 warning sites.** By instantiation count, the worst are the CORE cubed-sphere kernels:

| site | instantiations |
|---|---|
| `hydro/hydro_fluxes.cpp` 116, 308, 481 | 13 each |
| `mhd/mhd_fluxes.cpp` 122, 341, 568 | 8 each |
| `bvals/bvals_cc.cpp:111`, `bvals_fc.cpp:89`, `flux_correct_fc.cpp:201` | 1 each |
| `coordinates/coordinates.cpp` 475, 805, 901, 1005, 1036, 1092 | 1 each |
| `diffusion/viscosity.cpp` 120, 165, 204; `bvals_part.cpp` 358, 479; `track_prtcl.cpp:54` | 1 each |
| `bvals/flux_seam_cc.cpp:182` | 1 (mine, 806a9e6f) |

The gnomonic scheme is built on it: `pmy_pack->pcoord->GnomonicEquiangle{PrimFace,Flux,Emf}X*`
is CALLED from inside the flux kernels, and `pmy_pack->pmesh->GetPanelBoundary(...)` from
inside the halo pack kernels.

## Why it has never bitten

**Viper's GPU partitions are MI300A APUs with UNIFIED host/device memory**
([[viper-hip-build-recipe]]), so a host pointer really is addressable from a kernel. That
is why the dhj production GPU runs work. **On a discrete GPU (A100/H100/MI250X) this is an
illegal memory access**, so it is a portability bug, not a latent bug on this hardware.
There is no startup guard tying the cubed sphere to a particular device.

## Separately: cs_test.cpp does not compile for GPU AT ALL

Pre-existing, verified by compiling the file at **8671c912** with the same toolchain: three
static assertions from `size.template sync<HostMemSpace>()` / `mbpanel...` /
`ng...` at lines 1165, 1167, 1840 --
"Template parameter to .sync() must exactly match one of the DualView's device types or
one of the execution or memory spaces". So the cubed-sphere test pgen is **CPU-only**
regardless of the above, and every cubed-sphere result to date is a CPU result.

## How to check this yourself

The HIP configure in [[viper-hip-build-recipe]] is the only local way to compile-check
device code (login node `viper13` has no GPU, so it compiles but cannot run). 160 of 161
objects build; only `cs_test.cpp` fails. A brace-matching script that finds
`(pmy_pack|pmesh|pcoord|pmb|peos|pmbp)->` inside `KOKKOS_LAMBDA` bodies reproduces the
list without a compiler, but **prefer the compiler**: the script over-reports (nested
lambdas double-count, and pgen `if (pmbp->pmhd != nullptr)` hits inside host code look
identical).

## How it was fixed (af539941)

All the same shape: **give the kernel the values, not the object.**

* The eight gnomonic rotations left `Coordinates` for free functions in the new
  `coordinates/gnomonic_kernels.hpp`. They read just SIX Views between them, so
  `Coordinates::GnomonicTrigData()` bundles those into a `GnomonicTrig` POD the kernel
  captures by value. That alone killed the 13+8 instantiations in the flux kernels.
* `GetPanelBoundary` left `Mesh` for namespace scope -- it reads only its two arguments.
* Member Views and flags aliased into locals before the dispatch everywhere else.
* The reconstruction helpers became `static` and take their switches as arguments.
  `getWBq0`/`getWBerho` KEPT their member signatures as forwarders, so all 28 existing
  call sites (mostly pgens) were untouched -- worth copying when a static conversion
  threatens a wide blast radius.

**Two traps.** `static` CASCADES: a static caller cannot call a non-static callee, so
`GridPiecewiseLinearX1` dragged in `getWBq0` -> `getWBerho` -> `WBBackgroundStencil`, each
needing `wb_option` threaded through. And **never run `git stash` while a background build
is reading the tree** -- doing that during the lint parity check produced a phantom
"use of undeclared identifier GetPanelBoundary" that cost a debugging detour.
