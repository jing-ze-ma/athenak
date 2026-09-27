---
name: eventlog-mpi-deadlock
description: "FIXED 27da6380: <output>/file_type = log deadlocks ANY clean MPI run. WriteOutputFile's early return skipped the last_time advance, so 7 MPI_Allreduce ran every cycle and rank 0 blocked in a collective while neighbours waited on its sends. Likely affects upstream."
metadata: 
  node_type: memory
  type: project
  originSessionId: adaab265-737e-41c3-9d7d-e43d74682019
  modified: 2026-08-17T15:13:45.309Z
---

Found 2026-08-17 while chasing what looked like a general-EOS stall. It was not the EOS.

## The bug

`EventLogOutput::WriteOutputFile` began with

```cpp
if (header_written && no_output) return;
```

which returns BEFORE the `out_params.last_time += out_params.dt` at the bottom.
`last_time` is exactly what `Driver::Execute` (driver.cpp ~line 446) compares `time`
against to decide the output is due, so once that early return is taken the trigger stays
satisfied **forever** and `LoadOutputData` -- **seven collective `MPI_Allreduce` calls** --
runs on EVERY cycle thereafter.

It then deadlocks. Caught in gdb on four wedged ranks:

- rank 0: `PMPI_Allreduce` <- `EventLogOutput::LoadOutputData` <- `Driver::Execute`
- ranks 1-3: `MeshBoundaryValuesCC::RecvAndUnpackCC` <- `Driver::ExecuteTaskList`

Rank 0 reached the collective while its neighbours were still in the task list waiting on
its point-to-point sends, and a blocking collective does not progress those sends.

## Why nobody noticed, and why it looked like an EOS bug

**The counters are all zero exactly when nothing is going wrong.** So `no_output` is true
and the early return is taken on a CLEAN run. A run slamming into its density floor on 14%
of C2P calls writes a row every interval, advances the schedule, and sails through.

That is why it presented as general-EOS-specific: the general EOS does not hit its floors
on this problem, the ideal gas hits them constantly. It also cannot happen at the FIRST
write, which prints the header and so takes the normal path -- so it always looks like a
hang partway in.

**This is not general-EOS-specific and probably not fork-specific.** Any MPI run with
`<output>/file_type = log` that stops hitting floors will do it. Worth checking against
upstream IAS-Astrophysics.

## How it was diagnosed (the method is the reusable part)

1. `ndiag` was 2000, which hid everything. Setting `ndiag = 10` showed dt was NOT
   collapsing -- a flat 15.2288 -- so it was a hang, not a timestep death spiral.
2. Serial ran past the freeze point; every MPI rank count froze at the same cycle. That
   pointed at MPI, not physics.
3. A LOCAL mpirun on the login node reproduced it, which made `gdb -p <pid> -batch -ex "bt"`
   on each rank possible. Two backtraces settled it in one shot.
4. The freeze cycle (200) was exactly the second history-output time -- the giveaway that
   output scheduling was involved.

## Verified

4 ranks, general EOS, dhj 64x56x128: 500 cycles at a flat 0.665 s/cycle straight through
the old freeze point, mass and tot-E conserved, hst advancing on schedule.

## The other fix found on the way (`2a154e7d`) -- 21% faster

Not the stall, but real. `EOS_Data` has two spellings of several accessors; the 2-argument
ones SOLVE for the temperature (`SpecificHeatCv(d,e)` is documented "setup-time use only").
They were called per cell per stage with the temperature already in hand:

- `deep_hot_jupiter_rt`'s radiative source called `PresFromEint` then `TempKelvin` --
  two INDEPENDENT cold-start inversions of the same state -- then a Newton loop of up to
  100 iterations doing two more per iteration. Up to ~200 inversions per cell per stage.
  New `PresTempFromEint` gets both from one warm-started solve.
- `srcterms.cpp:181` / `srcterms_newdt.cpp:80` read `wtemp` two lines above, then asked
  `MeanMolecularWeight(d,e)` to solve for it again. **Generic** -- any general-EOS run with
  cooling.

**MEASURED 0.131 vs 0.165 s/cycle, 21%**, identical dedicated nodes.

## The rest of the audit, finished the same day (`d1c289dd`)

`wb_background.hpp` needed no external plumbing after all: every state in a stencil is
within a half cell of the anchor, so the anchor's temperature seeds all the others. Ten
inversions per stencil became ONE cold plus nine warm. Plus two refinements -- the outer
segments start from the interface temperature `WBAdvance` just returned, and
`WBEnergyFromEnthalpy` carries each Newton iterate's temperature into the next. The dhj
outer BC and `solar_convection`'s three SourceFunc sites read `wtemp`.

**Cumulative, MEASURED on the production config** (eos = general, general_eos = table,
ohmic_resistivity = eos, bbot = 10 G, dhj 64x56x128, 8x4, dedicated node, steady state
cycles 60-140):

| | s/cycle | |
|---|---|---|
| original | 0.1755 | |
| + RT/srcterms (`2a154e7d`) | 0.1406 | -19.9% |
| + wb_background (`d1c289dd`) | 0.1377 | -2.1% |
| **cumulative** | | **1.274x, -21.5%** |

The 2.1% UNDERSTATES the wb change: dhj reaches the well-balanced path only through its
outer boundary, ~3% of cells. `solar_convection` runs WBBackgroundStencil per cell inside
reconstruction, and that case is NOT measured.

**Verification, and a wrong hypothesis worth recording.** The ideal path is BITWISE
unchanged (history and every dump), as it must be since `Temperature()` ignores a guess
under a gamma law. The general path is NOT bitwise and cannot be: warm starting reaches the
same root by a different path and `SolveLog` stops on `|dz| < logtol`, so the converged T
moves by ~1e-13. dens/eint/bcc1/bcc2 are EXACTLY unchanged over 60 cycles; only
near-cancelling residuals move (velocities 5e-8 of rms, bcc3 1e-8). That signature is
identical to `2dedcbdb`, where the cause was FMA fusion and `-ffp-contract=off` made
everything byte-identical -- **so I ran that test here and it did NOT collapse the
difference.** Same signature, different cause. Do not assume.

FALSE POSITIVES, do not "fix": `srcterms.cpp:218` and `srcterms_newdt.cpp:116` are in the
relativistic cooling term, which the general EOS never reaches. One cold start remains in
`wb_background.hpp` by design -- the stencil anchor, which cannot be warm started without
threading `wtemp` through the reconstruction API.

Related: [[general-eos-table-cost]], [[general-eos-stage4-rho-e-table]],
[[dhj-ideal-xe-floor-relaxation]].
