# NOTE 2026-11-01 Caltech -> viper (USER): integ gates. (a) PASS on e32a104e and (a2) PASS on e1ecd210; (b) defaults block correct, BUT a fresh start spikes E_B 450x in cycle 1 when xe_tracer + ambipolar = eos + roche_wb_tan are all on. Do not launch r1m fresh on these defaults yet.

Answers TASK 13d5b766 and the e1ecd210 follow-up (NOTE 44cf9fd1).

## Binaries (H200, fmad off, incremental builds, -c 8)
| commit | binary | md5 |
|---|---|---|
| fin51-r-fix 8cf6db62 (reference for a) | athena_gpu_dhj_fin51_8cf6db62_fmadoff | (see .md5) |
| integ-test-1101 e32a104e | athena_gpu_dhj_integ_e32a104e_fmadoff | 4d9d2ecb19c20ba357e13d0f02fe4ba0 |
| d3f878be (reference for a2) | athena_gpu_dhj_d3f878be_fmadoff | 1b99f14bc69278203c407be29b1fca3c |
| integ e1ecd210 | athena_gpu_dhj_integ_e1ecd210_fmadoff | ba6abd1ce21cee1acd7d2557c06d55bc |

All runs: 1 rank, 1 H200 (hpc-sm-02-03 GPU 0/2), BVALS_CC_SPLIT=1.

## (a) Key-free bitwise, e32a104e vs 8cf6db62: PASS
keep/off_0p50.rst (developed 10 G r1m, cycle 8774), r1m keys from the rst, 150 cycles, hst every cycle.
- rst payload after `<par_end>`: IDENTICAL.
- dhj.mhd.hst: identical (150 rows).
- dhj.user.hst: only Efloor/Mfloor/Efloor_rt differ, <= 2.3e-13.

## (a2) Key-free bitwise, e1ecd210 vs d3f878be: PASS
Same setup. The rst payload is IDENTICAL and dhj.mhd.hst is identical. In dhj.user.hst (24 columns: Tdrag is added) only the floor columns differ, <= 1.4e-13.

## (b) Fresh r1m start with no special keys, e32a104e, 0.1 rot
- Deck: gates_1101/integ_b_sts.athinput, which is r1m_trc minus xe_tracer, floor_tracer, valfven_noinfall_fac, rkg_lean, rkg_dt_safety, top_bc and cs_seam_flux. ambipolar, roche_wb_tan and cs_seam_fc_nolimit were never in it. It keeps use_rkg_sts = true, ohmic eos, max_eta 1e15, valfven_max_cgs 2e6, roche true, am_hst_split false and flux_hst_floor.
- Start-up defaults block (e32a104e):
```
--- defaults (amdef-1019; <job>/recipe sets the default set: prod1019 | none) ---
  <mesh>/cs_seam_fc_nolimit = true (default on a fresh cs MHD start)
  <mhd>/am_transport = true (default on the cubed sphere; set it false to switch the angular-momentum transport correction off)
  <mhd>/valfven_noinfall_fac = 4 (default with valfven_max_cgs > 0)
  <mhd>/floor_tracer = true (default with valfven_max_cgs > 0)
  <mhd>/xe_tracer = true (default with ohmic_resistivity or ambipolar = eos)
  <mhd>/ambipolar = eos (default with ohmic_resistivity = eos and xe_tracer)
  <mhd>/rkg_lean = true (default on a uniform cubed sphere)
  <mhd>/rkg_dt_safety = 0.5 (default with STS on cs/sp)
  problem/top_bc = outflow1 (default)
  problem/roche_wb_tan = true (default with problem/roche)
  <time>/dt_retry = true (default)
--- end defaults ---
```
The same block on e1ecd210 (50-cycle fresh run) adds `problem/relax_drag = false (default on a fresh dhj start)`:
```
--- defaults (amdef-1019; <job>/recipe sets the default set: prod1019 | none) ---
  <mesh>/cs_seam_fc_nolimit = true (default on a fresh cs MHD start)
  <mhd>/am_transport = true (default on the cubed sphere; set it false to switch the angular-momentum transport correction off)
  <mhd>/valfven_noinfall_fac = 4 (default with valfven_max_cgs > 0)
  <mhd>/floor_tracer = true (default with valfven_max_cgs > 0)
  <mhd>/xe_tracer = true (default with ohmic_resistivity or ambipolar = eos)
  <mhd>/ambipolar = eos (default with ohmic_resistivity = eos and xe_tracer)
  <mhd>/rkg_lean = true (default on a uniform cubed sphere)
  <mhd>/rkg_dt_safety = 0.5 (default with STS on cs/sp)
  problem/relax_drag = false (default on a fresh dhj start)
  problem/top_bc = outflow1 (default)
  problem/roche_wb_tan = true (default with problem/roche)
  <time>/dt_retry = true (default)
--- end defaults ---
```
- Every requested default is listed. cs_seam_flux is absent from the block because "shared" is already its plain default since csedge-1019; it is not a recipe key.
- **Optional explicit-keys check: PASS.** No-key vs all keys explicit (integ_b_explicit.athinput), 100 cycles: rst payload IDENTICAL, mhd.hst identical, user.hst floor columns <= 7.9e-14.
- The 0.1-rot run itself finished clean: no FATAL, nan or NON-CONVERGED; 1757 cycles; mean dt 7.2 s (min 4.6); RKG 3.4 stages; 0.081 s/cycle.

### BUT: E_B spikes 450x in the first cycle
Total E_B from dhj.mhd.hst (cols 11-13):

| run | t = 0 | cycle 1 (t = 10.8 s) | 0.01 rot | 0.03 rot | 0.09 rot |
|---|---|---|---|---|---|
| (b) e32a104e, no keys | 1.797e31 | **8.23e33** | 2.19e32 | 1.86e32 | 1.55e32 |
| tg_fl (8ccae663, same deck with xe off and explicit production keys) | 1.797e31 | 1.80e31 | 1.80e31 | 1.85e31 | 2.56e31 |
| e1ecd210, no keys | 1.797e31 | **8.19e33** | | | |

- At 0.1 rot the (b) run has E_B 1.5e32 vs 2.6e31, KE 3.0e32 vs 6.9e32, and eos_dfloor 1.8x higher. That is a different start, not noise.
- 3-cycle fresh tests on e32a104e: the explicit-keys deck with ONE key changed on the command line, E_B at cycle 1:

| change | E_B at cycle 1 | spike? |
|---|---|---|
| none (all defaults) | 8.23e33 | YES |
| ambipolar = none | 1.88e31 | no |
| xe_tracer = false (so AD uses Saha x_e) | 1.797e31 | no |
| roche_wb_tan = false | 1.797e31 | no |
| cs_seam_fc_nolimit = false | 8.23e33 | YES |
| rkg_lean = false | 8.23e33 | YES |

- So the spike needs xe_tracer AND ambipolar = eos AND roche_wb_tan, all on together. Our T1 xe arms (2b72d789/89213fcb/ac98652f) had no AD and no wb_tan, and showed no spike.
- Cycle 0 log: `roche_wb_tan: equilibrium residual probe done at cycle 0` and `RKG super-stepping total stage=0 dt_diff=3.41` while dt = 10.8 s. The diffusive dt with AD is 3.4 s; it is 11.4 s without AD. Yet no RKG stages ran on cycle 0.
- My guess: on cycle 0 the wb_tan probe runs before the x_e tracer (or the AD dt) is set, so AD acts with an x_e ~ 0 eta_A (at the cap) without super-stepping. You know that code; please check.
- Restarts are unaffected: gates (a) and (a2) are bitwise, and the defaults only apply on fresh starts.
- **Recommendation:** fix before the r1m fresh launch. Until then, a fresh launch on these defaults needs `mhd/ambipolar = none` (or roche_wb_tan = false).

Run dirs (/resnick/groups/carnegie_poc/jingze/trc_1101/runs): ia_*, ib_sts, ib_nokey100, ib_expl100, isp_* (the 3-cycle tests), i2a_*, i2b_nokey50. Decks: gates_1101/integ_b_{sts,explicit}.athinput.
