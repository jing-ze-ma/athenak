---
name: fofc-compatibility-plan
description: "FOFC COMPATIBILITY PLAN 09-11 09:00 (user request: make fofc work with general EOS -> WB -> sp/cs -> MHD). No run uses fofc; no test covers it. Map of the 14 code paths and the ordered plan; decisions D1-D4 pending the user"
metadata:
  type: project
---
Findings (HEAD a10e367d): 1D segfault = hydro_fluxes.cpp:152-159 fofc widening lacks the one_d
branch (MHD has it) -> loops over k=-1..1 on 1-cell arrays. FOFC is called at the END of Fluxes
(hydro_tasks.cpp:426) after AddHeatFluxes/Viscous/AddGravFlux/RemoveWbFlux, and inside the WB
window where w0 holds the PERTURBATION (RemoveWbVar :295 .. AddWbVar :485) -> on WB runs the
fallback LLF sees rho-rho_bg (NaN) and the replaced flux re-adds the background pressure flux;
diffusive/gravity fluxes are deleted on flagged faces. Trial update hydro_fofc.cpp:61-79 uses
Cartesian dt/dx (wrong on sp; cs fatal-guarded at mesh.cpp:413; NO sp guard). sp_favg widening
is `else if` with fofc (hydro_fluxes.cpp:466,665) -> stale x2/x3 fluxes on sp+fofc. The c2p
only_testfloors pass tests dfloor||efloor||tfloor only: vceil never flags; on cs the deferred
floors (GnomonicEquiangleRaiseVel) are never run on utest. Scalar fluxes computed before FOFC
(hydro_fluxes.cpp:449) stay inconsistent with the replaced mass flux (every grid). MHD: no
transverse widening (mhd_fluxes.cpp:156), LLF fallback vs HLLE in cs_lowbeta_fallback, EMFs
overwritten; mhd_fofc resets flags by deep_copy.
Plan: step0 one_d branch + test_nr_fofc_sod_cpu; step1 general EOS (vceil flags FOFC, warm-start
guard, wder assert) + test_nr_fofc_geneos_cpu; step2 WB (full state in FOFC, re-apply WB/grav/
diffusive corrections on replaced faces) + test_nr_fofc_wb_cpu; step3 curvilinear (area/volume
trial update, additive sp widening, cs rotation + deferred floors, sp polar guard) + test;
step4 MHD (widening, fofc_rsolver option llf|hlle, EMF/CT check, vceil) + test_nr_fofc_mhd_cpu.
Invariants: fofc=false bit-identical to HEAD; fofc=true with no flagged cell bit-identical to off.
Decisions: D1 factor the cs deferred-floor block into an inline callable from the test pass
(recommended); D2/D3 move FOFC BEFORE all additive flux corrections (fallback = pure hydro flux
of the FULL state, then WB removal/gravity/diffusion apply to it too) (recommended); D4 fix the
scalar-flux inconsistency (breaks bitwise only for nscalars>0 && fofc).

## DECIDED 09-11 09:20 (user "ok"): call order = FOFC right after the Riemann fluxes, before every
## additive correction (published AthenaK semantics: trial update from Riemann fluxes only, no
## sources; fallback DC + LLF/HLLE; before AMR flux correction -- arXiv 2409.16053, 2409.10384,
## 2609.06150); WB fallback from the FULL state; scalars recomputed from the new mass flux; cs
## deferred floors factored out for the test pass. Agent doing steps 0+1+call order + tests
## test_hydro_fofc_sod_cpu / test_hydro_fofc_blast_cpu; invariants I1 (red giant GV3b) + I3 (GR).

## STEPS 0+1+CALL ORDER DONE, committed 09-11 10:20 (local, not pushed): one_d branch; FOFC before
## all additive corrections (hydro+MHD); vceil_test flags on the test pass; general-EOS test pass
## bails on non-finite; wder fatal; WB full state via wfull accessor; scalar fluxes recomputed;
## tests tst/test_suite/hydro/test_hydro_fofc_{sod,blast}_cpu.py (7 pass; inputs in tst/inputs/).
## Invariants: red giant GV3b bit-identical, GR monopole bit-identical, HEAD segfaults on 1D fofc.
## LEFT: MHD has NO vceil at all; mhd_fofc has no WB full state; WB test not written; step 3
## curvilinear (area/volume trial update, additive sp widening, cs rotation + factored deferred
## floors, sp polar guard); step 4 MHD. Regress builds b_edit/b_head/b_fofc_rg.

## STEP 3 CURVILINEAR DONE 09-11 (UNCOMMITTED in tree, verified by me): hydro_fofc.cpp area/volume
## trial update (no geometric sources, like RKUpdate), cs rotation of the single-state solves via
## ldx*/rotflx* (mirrors gnomonic_kernels PrimFaceX*/FluxX*), FOFC-csfloor kernel runs the factored
## deferred floors (new coordinates/gnomonic_raisevel.hpp, GnomonicEquiangleRaiseVel now calls it) on
## the trial state; mesh.cpp cs fatal now mhd/fofc only; sp_favg widening = union (was a no-op). Polar
## axis safe by area.x2f=0; Debug bounds-check clean. Tests tst/test_suite/hydro/test_hydro_fofc_{sp,cs}
## _cpu.py + tst/inputs/fofc_{sp,cs}.athinput (9 fofc tests pass). I1: sp/cs ideal, cs general,
## cs general+vceil, Cartesian fofc=true all DATA IDENTICAL vs HEAD worktree
## /orion/u/jinma/ATHENAK/wt_fofc_head (remove when done). Builds build_fofc3, build_fofc3_dbg.
## NOTE: coordinates.cpp now mixes step-3 refactor with the uncommitted red-giant
## dfloor_keep_temperature hunk (moved verbatim into the header) -> split carefully at commit.
## --style fails only on ~2200 pre-existing violations. LEFT: step 4 MHD.

## STEP 4 MHD DONE 09-11 (UNCOMMITTED, user said HOLD COMMITTING): mhd_fofc.cpp sp area/volume trial
## update + spherical curl for the trial bcc (old form missed r and sin theta) + wfull WB accessor
## (MHD has RemoveWbVar too); <mhd>/vceil (shared EOS_Data::vceil; ideal+general MHD c2p; refused
## for SR/GR and cs); <mhd>/fofc_rsolver=llf|hlle (hlle_mhd_singlestate); widening was already
## sufficient (comment only); EMF/CT consistent (FOFC overwrites face EMFs before CornerE).
## Test tst/test_suite/mhd/test_mhd_fofc_blast_cpu.py (21 pass), inputs fofc_mhd_blast/fofc_mhd_sp.
## I1 fofc=false Cartesian+sp identical vs HEAD; I2 fofc=true unflagged identical. Build build_fofc4.
## Harness needs gcc 13.1.0 CXX + anaconda python>=3.7. LEFT: mhd/fofc on the cubed sphere (fatal
## kept: needs gnomonic rotation of single-state solves, trial curl, EMFs; vceil via RaiseVelMHD).
## ALL FOUR STEPS DONE; the whole FOFC stack is uncommitted on top of the red-giant RT hunks.

## STEP 5 MHD ON CS DONE 09-11 (UNCOMMITTED, HOLD COMMITTING): mhd_fofc.cpp ldx*/rotflx* with the
## field (FaceBX2 on x2 loads, FluxX* rotate-back, EmfX1 on the replaced e3x1 -- the only EMF the
## high-order sweep rotates, verified), cs trial curl (approximate: edge lengths, documented),
## new coordinates/gnomonic_raisevel_mhd.hpp factored from GnomonicEquiangleRaiseVelMHD (now a
## parallel_reduce counting neos_vceil) with <mhd>/vceil; cs refusal lifted in eos.cpp; mesh.cpp
## fatal gone. shock_tube.cpp gained <problem>/bazi (curl-of-A azimuthal field on cs, inert at 0).
## Test tst/test_suite/mhd/test_mhd_fofc_cs_cpu.py (6 pass, 2.5 min) + inputs/fofc_mhd_cs.athinput;
## its floor-reduction assert is only ~8% (physical floor). I1 identical vs HEAD; I2 identical.
## Build build_fofc5. FOFC PLAN COMPLETE: 41 fofc tests pass. Uncommitted diff = steps 3+4+5 on top
## of the red-giant RT hunks (eos/*, hydro.cpp, two_stream_rt.hpp, coordinates.cpp keep_temperature).

## WB + MPI/SMR TESTS DONE 09-11 (UNCOMMITTED): new built-in pgen src/pgen/tests/wb_atm.cpp
## (registered in pgen.hpp/pgen.cpp/CMakeLists; no built-in pgen filled u0wb/w0wb before) +
## tst/inputs/fofc_wb.athinput + test_hydro_fofc_wb_cpu.py (4 pass): full-state fallback proven
## (WB-on vs WB-off L1 ratio 1.2; forcing the perturbation fallback -> NaN by cycle 44).
## test_hydro_fofc_blast_mpicpu.py + test_mhd_fofc_blast_mpicpu.py (2+2 pass): 4 ranks, 16 root
## blocks + static refinement, nghost=4; 1 vs 4 ranks bit-identical incl. FOFC counts.
## MPI harness: module load gcc/13 openmpi/4 on the same command line. WB scheme itself NaNs above
## ~Mach 8 kicks regardless of FOFC (known "less robust" note in hydro.hpp). 49 fofc tests total.

## DYNAMIC WB CASES ADDED 09-11 (UNCOMMITTED): test_hydro_fofc_wb_cpu.py parametrised
## static/dynamic (8 pass); wb_atm.cpp fills phicc0/phi0 faces for dynamic and uses the local
## background face-pressure difference as the momentum source (as hse_atm/red_giant). Dynamic
## balanced state holds to 1e-16, FOFC flags 0. MACH CEILING: static+reconst_perturb NaNs above
## ~Mach 8 regardless of FOFC; DYNAMIC WB fofc=false NaNs above Mach ~39, fofc=true survives to
## Mach 80 (every amplitude tried). 53 fofc tests total.
