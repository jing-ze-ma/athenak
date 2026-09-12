---
name: cs-gpu-multiblock-fault
description: FIXED 4302a008 -- the cubed-sphere multi-block GPU fault was CAPTURED STATE, not logic: one kernel named both sendbuf[56] and recvbuf[56]. Split by rank so each captures one. Records the three-way A/B that exonerated the cross-block write, and that an INERT change hid the bug completely
metadata:
  type: project
---

**FIXED 2026-08-29, 4302a008.** Found on `apudev` (MI300A); I introduced it in 806a9e6f.
Cubed sphere now runs on GPU at 1, 4 (24 blocks) and 16 (96 blocks) MeshBlocks per panel,
hydro and MHD, and on 2 GPUs under MPI with 1-vs-2-rank global sums BIT-IDENTICAL.

## The symptom

With **2x2 MeshBlocks per panel** (24 blocks), a cubed-sphere run on one GPU dies with

    Memory access fault by GPU node-2 (Agent handle: ...) on address 0x... Reason: Unknown.

**11 of 20 trials**, hydro and MHD alike. Interleaved A/B in a single job, so the
comparison is not confounded by node or session.

## What is established

| | |
|---|---|
| 1 MeshBlock per panel | **0 / 20** -- clean, and correct (matches CPU) |
| 2x2 or 4x4 per panel | faults intermittently |
| CPU, every decomposition | clean; gates at round-off, NULL control exactly 0 |
| Kokkos-debug / bounds-check build | **never faults**, and bounds checking NEVER fires |
| faulting dispatch | `MeshBoundaryValuesCC::PackAndSendFluxSeamCC`, via `AMD_SERIALIZE_KERNEL=3 AMD_LOG_LEVEL=4` |
| buffer indices | **audited on the host for 1/4/16 blocks per panel: no overrun** |

**THE METHODOLOGICAL TRAP, and the most reusable thing here.** The bug is latent enough
that a semantically INERT extra branch in `SeamFaceGeom` took it from 11/20 to **0/20**
by changing code layout alone. I nearly reported that as the fix. **Never accept "the
symptom stopped" on this code**: interleave the A/B in a single job, run 20 trials, and
require a mechanism that explains WHY.

## THE CAUSE: captured state, not logic

`sendbuf` and `recvbuf` are `MeshBoundaryBuffer[56]`; each element holds several Views
and six index triples, so each array is tens of kilobytes. The pack kernel NAMED BOTH, so
`[=]` copied both into one functor -- large enough that Kokkos launches it from GLOBAL
memory instead of constant memory (visible in the ShaderName), and evidently too large.
**Fix: two kernels, on-rank capturing `recvbuf` alone and MPI capturing `sendbuf` alone.**
The send and receive index ranges are identical for a seam face, so each kernel can read
them from whichever array it already holds.

**This is the SECOND bug of exactly this shape in one session** -- the first was handing
the flux kernels all six gnomonic trig Views ([[gpu-this-capture-device-lambda]]). The
rule: **give a kernel only what it uses**; on this code an over-large functor does not
fail cleanly, it faults at a wild address.

## The three-way A/B that found it (interleaved, one job, 20 trials each)

| variant | faults |
|---|---|
| both arrays captured, cross-block write ON | **12/20** |
| `sbuf` only, cross-block write OFF | 0/20 |
| `rbuf` only, **cross-block write ON** | 0/20 |

The third row EXONERATES the cross-block write, which was my prime suspect. Without it I
would have "fixed" the wrong thing.

## Hypotheses RULED OUT -- do not re-test these

* **Cube-vertex collision** (the FC path's old bug, [[cubed-sphere-seam-emf]] fix 3):
  refuted. The destination map is a clean bijection with full reciprocity -- 0 collisions,
  0 non-reciprocal -- at 1, 4 and 16 blocks per panel.
* **Buffer overrun**: refuted. A host audit of every write, all three decompositions.
* **`dest`/`dm` out of range**: refuted, same audit.
* Kokkos bounds checking never fires; a Kokkos-debug build never faults.

Harness: `/viper/u2/jinma/cs_gpu_check/` (outside the repo) -- `job*.sh`, `athena.clean`
(pre-classifier), `athena.sfg` (post), `athena.pre` (806a9e6f), `athena.bc`, `athena.dbg`.
`apudev` is 15 min, 1 node, 2 GPUs. **The scratchpad on the login node's /tmp is NOT
visible to compute nodes** -- submit from `/viper/u2/...`.
