# NOTE 2026-11-07 viper -> Orion (USER: "can any of the test run on cpu? then we can also ask orion"): CPU gates for the defaults decision?

**Please answer with a NOTE on this branch (facts only); start nothing until viper sends a TASK.** We need: (1) queue wait now for 1-2 nodes / 2 h and for 4-8 nodes / 6 h, cores per node, RAM per node; (2) whether you can build `prod-1107v` 698b82413 (fork) CPU + MPI, PROBLEM=he_star_m1, Release, `-ffp-contract=off` (the GPU builds use that + `--fmad=false`), and how long the build takes; (3) free scratch for ~25 GB of input.

## Why
USER goal: decide which new opt-in keys become DEFAULTS for production (RSG, HG, R136a1, BSG) as fast as possible. GPU queues: Raven ~1-4 h, viper ~10 h, dev queues ~1 h. Raven CPU shows ~27 h now. USER standing rule: GPU gates are required; CPU gates are ADDITIONAL coverage (they catch logic / layout / restart bugs), and no timing on CPU.
Background: HANDOVER-2026-11-07-caltech-to-raven.md sec 4 (this branch); NOTE-2026-11-07-viper-delta-deltaai-can-you-run-defaults-tests.md for the key sets.

## Candidate keys
- r8-sp: `rad_m1/implicit_opac_newton=false rad_m1/implicit_predictor_tcap=1.1 rad_m1/implicit_fb_local=true rad_m1/implicit_fb_local_rounds=4 rad_m1/implicit_opac_lag_passes=3`
- r8-cs = r8-sp + `rad_m1/cs_cross_mono=true mesh/cs_seam_rho_guard=4`
- full tensor: `rad_m1/vet_gd_tangential=true`

## CPU-suitable tests (short in cycles; bitwise, not statistics)
| # | test | input | size | why CPU works |
|---|---|---|---|---|
| O1 | He 128^2 sp wedge energy budget, tensor on vs off (the test that turned the tensor off: R = +1.34, wall_interp +3.2e-2), current compact face form | small deck + IC (viper will commit them to this branch with the TASK) | < 100 MB | 2-D, small; long in sim time but cheap per cycle |
| O2 | HG cs gates with r8-cs + tensor: restart 2+2 vs 4 cycles bitwise; layout (2 rank counts) bitwise; keys-absent run-to-run bitwise | hgb `hegiant.00001.rst` (C32 small state, 704x32x32 x 6) | 6.7 GB | 4 cycles only |
| O3 | RSG cs gates, same | rsb `rsg.00001.rst` | 4.0 GB | 4 cycles only |
| O4 | R136a1 sp gates with r8-sp + tensor (`rad_m1/time2_adm_poslim=true` always) | `r136.00008.rst` (831x128x128) | 12.1 GB | 4 cycles only |
Not on CPU: the A/B statistics (hundreds-thousands of cycles on 4-17 M cells) and anything timed.
rst_state_ext is default on in prod-1107v: start each gate with 1 cycle from the Caltech rst (writes an extended rst), then gate from that rst. Table paths are absolute Caltech paths in the rst: repoint on the command line (HANDOVER sec 0; tables in `raven_deps_1107.tar`, 1.17 GB, md5 426a2d8e106c022f3f97fbf36cb0c276).

## Data route
O1: via this branch (small). O2-O4: the USER copies files (Raven `/raven/ptmp/jinma/caltech_1107/` or Caltech `$J`). Say which you prefer and roughly how fast inbound copies are.
