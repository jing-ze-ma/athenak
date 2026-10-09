# NOTE 2026-10-09 (viper): default-on of the 10-09 production set on rt-integration (commit 40f003cf)

The user asked for this on 10-09 ("allow the twin default edit", "make ic_eint_from_t default"). Only effective
defaults change; every key is still read, and naming it wins.

## Rules for every new default
- The new defaults apply to FRESH runs only. On a restart whose file lacks a key, the old default is kept.
- Every resolved value is recorded in the parameter dump.

## New defaults

| key | new default | where it applies (elsewhere the old default) |
|---|---|---|
| rad_m1/vet_gd_twin | true | vet_gd runs, unless vet_gd_async = true |
| rad_m1/vet_gd_twin_fuse | true | with the twin, where VetGdInit's needs hold: vet_gd_iter 1, no async / band_exit, and with MPI vet_gd_halo_compact > 0; else false (the unfused twin) |
| rad_m1/vet_gd_twin_det | true | with the twin |
| rad_m1/vet_source_noesrc | true | with the twin; acts under time_scheme = be only |
| rad_m1/vet_gd_twin_lowmem | false (unchanged) | |
| rad_m1/implicit_face_opac_n | true | all implicit runs; acts on blend faces only |
| rad_m1/implicit_pos_floor_solve | = implicit_pos_floor (true on fresh runs) | |
| rad_m1/report_newton_fb | false (already the default) | |
| rad_m1/implicit_precond | mg also on the spherical-polar wedge | fast path, multi-D, single level, no polar boundary, no cubed sphere, no krylov_dev; non-fast path stays line; cs / polar / SMR stay rbgs_fwd |
| rad_m1/implicit_flux | blend | multi-D transport = implicit with a VET closure that carries ray data (closure = vet_sc, or vet_col + vet_gd = true), not cubed sphere |
| rad_m1/implicit_flux_faces | all | the same, and only when implicit_flux resolves to blend |
| rad_m1/implicit_flux_beam | halfrange | the same, and only when faces resolve to all |
| rad_m1/implicit_blend, _tau0, _flo, _fhi, _xthin | idort_f, 1, 0.3, 0.6, 30 | the same, and only when the beam resolves to halfrange |
| rad_m1/vet_scatter | true | VET with ray data (as above), opacity = table and <hydro\|mhd>/general_eos = table |
| rad_m1/vet_scatter_kappa_e | eos | wherever vet_scatter is on and general_eos = table; else const (form ma and j relaxed were already the defaults) |
| hydro/sp_x2_periodic_image (and mhd/) | true | spherical-polar, multi-D, theta-periodic, not theta-stretched (the only place it acts) |
| problem/he_ic_eint_from_t | true | he_star_m1 with he_ic_cols = 5 and an IC file with >= 7 columns on its first data line; else false |

Notes:
- M1, eddington, tau, and vet_col without vet_gd keep the old scheme defaults.
- An input that names implicit_flux = central keeps faces x1, beam closure and blend tau_f.
- he_ic_eint_from_t reads column 6 as T [K]. A data line without columns 6-7 is a FATAL (existing check), so the
  pgen never silently misreads a column.

## Gates
Raven gpudev 4 A100, job 31034693 (build 31034692, binary
/raven/ptmp/jinma/merge_1009/bin/athena_he_a100_40f003cf_dflt1009). The reference is the production binary
98835d99 on the production input; bin data were compared with /viper/ptmp2/jinma/merge_1009/bincmp.py.

| input | all explicit keys kept | keys equal to the new defaults removed |
|---|---|---|
| AG Car A hr, 10 cycles | hst 2/2 byte-identical, bin data 24/24 bitwise | the same, 24/24 |
| He giant fresh N897 eft, 3 cycles | the same, 15/15 | the same, 15/15 |
| BSG true-repro hr (reduced 4 blocks), 5 cycles | DIFFERS (18/28 dumps), as expected: he_ic_eint_from_t is now true for its 7-column IC | with the pin he_ic_eint_from_t = false: hst identical, bin data 28/28 bitwise |

The keys removed are listed by /raven/ptmp/jinma/merge_1009/mk_r.py; the generated inputs are in
/raven/ptmp/jinma/merge_1009/dgate.

CPU checks (login node; results in /viper/ptmp2/jinma/merge_1009/cpugate/RESULTS.txt):
- WASP-121b cubed-sphere MHD: byte-identical to the merge build, all 8 files (hst and bin).
- He box M1 test (vet_sc, central flux, 4 ranks): runs with no FATAL and CHANGES, as intended. On it vet_scatter
  turns on with kappa_e = eos (table opacity + general EOS), and implicit_pos_floor_solve and implicit_face_opac_n
  turn on.
- Regression subset (rad_m1, hydro, mhd, rad; regress_dflt2.log):
  - `_cpu`: 62 passed, 3 skipped.
  - `_mpicpu`: 4 passed, 3 skipped, 1 failed. The failure is rad/test_rad_lwave2d_amr_mpicpu, which segfaults on the
    old code too.
  - No test threshold moved, so no pins and no threshold changes were needed in tst/.

## Pins: how to reproduce pre-10-09 behaviour
Name these in the input; the old values are shown.

- BSG true-repro (DONE: pinned, comment "PIN (defaults-1009)"): `<problem> he_ic_eint_from_t = false`.
  - viper /viper/ptmp2/jinma/bsg_truerepro_1008/run3d/bsg3d_truerepro2.athinput and _hr.athinput.
  - Raven /raven/ptmp/jinma/truerepro2_1009/inputs/bsg3d_truerepro2_raven.athinput and _raven_hr.athinput.
  - Templates on fork/bsg-files-1009 ddf31ea6: bsg3d_truerepro2.athinput.in, bsg3d_truerepro2_hr_lm.athinput.in,
    bsg_hr_dc5.athinput.in.
  - Caltech's BSG copy runs binary 6c5d8fb2 and is unaffected.
- AG Car A/B and He giant: no pin needed. They already name every new default.
- He box M1 / any vet_sc or vet_gd input that does not name them: name these keys with their old values.
  - `<rad_m1>`: vet_scatter = false, implicit_face_opac_n = false, implicit_pos_floor_solve = false.
  - vet_gd only: vet_gd_twin = false (the twin keys and vet_source_noesrc then stay off).
  - Under vet_sc / vet_gd without a named implicit_flux: implicit_flux = central.
  - sp wedges: implicit_precond = rbgs_fwd (sp fast path).
  - Theta-periodic sp: `<hydro>` (or `<mhd>`) sp_x2_periodic_image = false.
  - he_star_m1 with a 7-column IC: `<problem>` he_ic_eint_from_t = false.
- WASP-121b / dhj cubed sphere: nothing changes and no pin is needed (no rad_m1, not spherical-polar). Gate:
  byte-identical, above.
- Restarts: unaffected unless the restart file lacks the key AND the run is fresh. Restarts never take a new default.
