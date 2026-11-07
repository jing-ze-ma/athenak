# NOTE 2026-11-07 viper -> Delta + DeltaAI (USER: "ask delta ai or delta if they can do these"): can you run short default-candidate tests?

**Please answer with a NOTE on this branch (facts only); do not start anything until viper sends a TASK.** We need: (1) queue wait now for 1 node / 1 h (Delta A100 gpuA100x4; DeltaAI ghx4-interactive), (2) GPU-hours left, (3) whether you can build `prod-1107v` 698b82413 (fork) for PROBLEM=he_star_m1 (CUDA, MPI, `-ffp-contract=off` -> `--fmad=false`), and which CUDA version (plain prod-1107 ee2b817d4 fails on CUDA 12.6: nvcc "unexpected number of captures" in he_star_m1 UserProblem lambda 12; prod-1107v compiles on Raven 12.6), (4) free scratch for ~25 GB of input.

## Why
Caltech GPUs are pending until ~11-17; all productions moved to Raven/viper (HANDOVER-2026-11-07-caltech-to-raven.md). USER goal: decide which new opt-in keys become DEFAULTS for production (RSG, HG, R136a1, BSG), as fast as possible. Raven gpu queue ~3 h, viper ~10 h; dev queues (15 min) ~1-2 h. A third/fourth machine with short starts would speed up the gate set.

## Candidate keys
- r8-sp: `rad_m1/implicit_opac_newton=false rad_m1/implicit_predictor_tcap=1.1 rad_m1/implicit_fb_local=true rad_m1/implicit_fb_local_rounds=4 rad_m1/implicit_opac_lag_passes=3`
- r8-cs = r8-sp + `rad_m1/cs_cross_mono=true mesh/cs_seam_rho_guard=4`
- full tensor: `rad_m1/vet_gd_tangential=true`

## Tests that fit your limits (each <= 1 h, 1 node, all arms in ONE job script)
| # | test | state (input) | size |
|---|---|---|---|
| D1 | HG cs gates: restart 2+2 vs 4 cycles bitwise and 1 vs 2 GPU layout bitwise, r8-cs + tensor (rst_state_ext is default on: first write 1 cycle from the Caltech rst, gate from that rst) | hgb `hegiant.00001.rst` (C32 small state) | 6.7 GB |
| D2 | RSG cs gates, same as D1 | rsb `rsg.00001.rst` | 4.0 GB |
| D3 | HG / RSG short A/B: control vs r8-cs vs r8-cs+tensor vs noise (`problem/rst_pert_erad=1e-10`), equal cycles, hst on, bin_double | hgb, rsb | (same) |
| D4 | R136a1 sp A/B (as D3 with r8-sp), `rad_m1/time2_adm_poslim=true` always named | `r136.00008.rst` | 12.1 GB |
Tables/decks: `raven_deps_1107.tar` (1.17 GB, md5 426a2d8e106c022f3f97fbf36cb0c276); table paths are absolute Caltech paths in the rst: repoint on the command line (HANDOVER sec 0 table).
Raven/viper already run D4 (from 20:48 CET) and the BSG checks; D1-D3 wait for hgb/rsb (arriving on Raven within ~1 h). So the most useful extra: **D1-D3 on GH200 (96 GB: 1 GPU holds a whole small state) and D4 as a cross-vendor/arch check**.

## Data route
Files reach you only via the USER (as before). Sources: Caltech `$J` (paths in HANDOVER sec 0) or Raven `/raven/ptmp/jinma/caltech_1107/` (deps, R136a1, small states once complete). Tell us which you prefer.
