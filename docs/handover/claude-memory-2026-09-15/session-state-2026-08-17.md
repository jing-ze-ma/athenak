---
name: session-state-2026-08-17
description: "START HERE. End of 2026-08-17: origin/general-eos = d1c289dd, four commits pushed. Floors/max_eta/STS settled and documented, an MPI deadlock fixed, general EOS now 3.07x ideal not 4.4x. Stage 4 re-baselined and still not started."
metadata: 
  node_type: memory
  type: project
  originSessionId: adaab265-737e-41c3-9d7d-e43d74682019
  modified: 2026-08-17T15:48:37.605Z
---

End of 2026-08-17. Supersedes [[session-state-2026-08-16]] as the entry point.

## ALL PUSHED -- `origin/general-eos` is at `d1c289dd`

Working tree clean under `src/ inputs/ tst/ docs/`. Four commits today, oldest first:

- `6003c1ce` `docs/ideal_gas_resistive.md` + `inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput`
  -- the recommended settings for ideal-gas resistive runs ([[dhj-ideal-xe-floor-relaxation]])
- `27da6380` **the event log deadlocked every clean MPI run** ([[eventlog-mpi-deadlock]])
- `2a154e7d` stop re-solving for a temperature the code already has -- 19.9%
- `d1c289dd` the same for the well-balanced background -- 2.1% more

Also pulled at the start of the day: `d52f3610` (GPU host-pointer fix) and `0accdbdb`
(meshblock porting advice), both authored on Viper.

## What changed in the numbers everyone quotes

**The general EOS is no longer 4.4x.** Measured with `ohmic_resistivity = eos` on both
sides, identical input and hardware, restricted to the window where BOTH are hydro-limited:
**3.07x per cycle, 2.32x per simulated second.** Details and the trap in
[[general-eos-table-cost]].

**Recommended ideal-gas resistive settings changed too:** `max_eta = 1e13` with
`use_rkg_sts = false`, 1.54x faster than 1e14+STS at the same dt. STS is NOT a standing win
at 1e14 -- the sign flips mid-run as dt_hydro decays.

## The three traps found today, all recorded in their own notes

1. **Never diagnose a hang with `ndiag = 2000`.** It hid that dt was healthy and turned a
   deadlock into a mystery for hours. See [[eventlog-mpi-deadlock]] for the method that
   worked: serial vs MPI, then a LOCAL mpirun so `gdb -p` was possible.
2. **`--cpus-per-task` must stay on p.exclusive** whenever OMP threads are used. Dropping
   it (as [[freya-job-submission]] used to advise) oversubscribes and looks exactly like a
   code pathology -- it read as a "140x slowdown in the general EOS" and was nothing.
3. **The same numerical signature can have different causes.** dens/eint/bcc1/bcc2
   identical with velocities and bcc3 moving at ~1e-8 was FMA fusion in `2dedcbdb` and is
   NOT in `d1c289dd`; `-ffp-contract=off` distinguishes them and was run both times.

## Thread status

| thread | state |
|---|---|
| ideal-gas resistive: floors, max_eta, STS | **CLOSED**, shipped ([[dhj-ideal-xe-floor-relaxation]]) |
| event log MPI deadlock | **CLOSED** ([[eventlog-mpi-deadlock]]) -- but check upstream, it is not fork-specific |
| cold-start EOS call audit | **CLOSED**, one deliberate anchor solve left |
| general EOS Stage 4 (rho,e) table | **RE-BASELINED, NOT STARTED** -- worth ~1.5x now, not ~2.1x ([[general-eos-stage4-rho-e-table]]) |
| solar_convection atmospheric runaway | **STILL OPEN**, untouched today ([[solar-convection-general-eos]]) |

## LEFT RUNNING AT SESSION END -- do not cancel

Four long dhj runs in `idflr/L_{id13,id14,gen13,gen14}`: ideal vs general x max_eta
1e13 vs 1e14, all WITHOUT RKG, whole node each, tlim = 1e7, 24 h wall. Two were running
healthily and two still queued. Full details and the analysis recipe in
[[dhj-ideal-vs-general-cost]] -- **read that first next session.**

## Immediate next steps if picking up Stage 4

The user asked to start it, then asked for the stall first; the stall is now fixed, so
Stage 4 is unblocked and the benchmark configuration actually runs. Task 1 is unchanged and
still not done: **compute the real e range over the intended (rho,T) box with a host-side
harness calling `EOSCompositionModel`**, to see how non-rectangular the (rho,e) domain is.
A throwaway Saha model is NOT good enough -- one was tried and was visibly wrong at the cold
end.

Also unmeasured and worth a run: `solar_convection` with `d1c289dd`, where
WBBackgroundStencil runs per cell in reconstruction rather than only on ghost cells, so the
2.1% measured on dhj should be much larger.

## Scratch directories on /orion (none in [[run-directory-untouchable]])

- `idflr/` -- everything from today. `base d3 d6 all3 p3 nosts eta13` (the floor and
  max_eta/STS sweep), `q_ideal q_gen` (the EOS cost comparison), `p_orig p_before p_after`
  (the fix attribution), `evtest` (the deadlock reproducer), plus binaries
  `athena_ref` (pre-fix), `athena_evfix`, `athena_wb` (current), `athena_proxy`
  (`logtol = 1e30`, the no-inversion proxy) and `athena_fpc_{old,new}` (-ffp-contract=off).
- `build_solar/` (PROBLEM=solar_convection), `build_fpc/` (-ffp-contract=off),
  `athenak/build/` -- now PROBLEM=**deep_hot_jupiter_rt**, no longer solar_convection.
- Older: `soleos/ dhjb/ xetest/ build_tests/ build_dhjrt/`.
