# NOTE 2026-11-24zzzr orion -> Caltech: HG gate for rad_vet/time2_gtime_slope (USER approved; run after the 24zzzh batch and the 00205 copy)

Branch pushed: `caltech-sdorder-a840d7e92` (TEST branch, not for merge). Args: docs/handover/hg_gate/IA123_hg_restart.args (MERGE123; checked only on a small old-scheme HG rst). hegiant.00205 from viper per 24zzzk (check its md5). HG on T2 has never run on H200: do the +1 step check first and post peak GPU memory. Run when the test batch leaves 2+ nodes free; nothing on production. CPU state: the key restores order 2 and the shock errors on Orion CPU; HG is the open risk (see FAIL rules).


**Why.** In the release, sdirk2's exchange is first order in time where it is stiff (wave p 1.0-1.5, M5 shock rho L1 +27 % over
hesdirk2). The key restores order 2 in the CPU gates. The risk is HG: the old f = 1/g (HEUNFIX FHD2) dropped HG's dt from 90 s to
27-51 s after ~14 cycles. This key never adds more than that form, and adds nothing where the exchange changed sign, but a spike
that keeps its sign passes. On the CPU Marshak front the flux state's max |ef/es - 1| grew from 6e-3 (off) to 0.8 (on).
HG is the gate that decides.

## Binary
One TEST binary for every arm: branch vet-sdorder-test a840d7e92 = release 89b05b1c7 + 190e1b88e (the key) + the SDO_PRINT diagnostic
(diff: stag2/sdorder/sdo_print.diff). PROBLEM = he_star_vet, your CUDA 12.9 + g++ 13 recipe, .prov checked. With the key off and
SDO_PRINT unset it is the release (Orion CPU: key-off data-bitwise, 258 files of the quick suite; NOT checked on H200).

## Start, keys, layout
- rst hegiant.00205 (md5 82058331..., viper/Raven copy now at Caltech; N0 794784). 24 MeshBlocks 704 x 64 x 64: 4, 6, 8 or 12 ranks.
  The same rank count in every arm. Memory unknown at 4 H200: use 8 ranks if 4 does not fit.
- Keys = IA123_hg_restart.args + `rad_vet/vet_gd_allow_clamp=false` + `SDO_PRINT=1` in the environment of EVERY arm
  + `job/final_restart=false` + hst every cycle (output1/dt=0 output1/dcycle=1), no bins.
  allow_clamp: not stored, not named in the recorded IA9 HG gate (resolved/hgr_IA9_HGRk.input.txt), so it ran the default false = FATAL on a clamped read.
- Compared with the recorded IA9 gate (resolved/hgr_IA9_HGRk.input.txt): conv_max 1e-5, change_form linear, xdonor_freeze 4, fofc true and
  rad_signal_speed true are the SAME (freeze, fofc and signal speed from defaults / stored keys; hg_latest's 1e-4 is an older set).
  DIFFERENT from the record, so key off is NOT the record: rad_signal_speed_diff / _phase = true (this changes dt), fofc_retest,
  floor_pay_skip, dfloor_flux_limit, floor_from_ekin, layout tile 64, clean_stag_central, clean_stag_eicarry, clean_pred_xch,
  clean_clos1_budget, clean_bud_maxnorm, clean_gas_secant = max_stage2, fwb_stage / fwb_t2 = false, w_fse_tau, wb_exch_local,
  vet_rotate_avg_ema_wb = true, clean_stag_mg_levels 2, gtime_mom, rotate_every / t0 / stretch. FHD2C / FHE0 used the FF keys of the Raven
  hgrep harness, which I could not read (Raven down), so their conv_max and change_form are NOT established here.
  Therefore score the PAIR; the record is orientation only.
- **One check before the arms:** a +1 step with `time/nlim=794785`: rc 0, no FATAL, CLEAN WARNING 0. If it FATALs (ema_wb vs the file's stored
  value is one candidate; MERGE123 checked the args only on the small HG rst), stop and report.

## Arms (one at a time or on separate nodes; nothing on production nodes)
| arm | key | cfl | cycles |
|---|---|---|---|
| H0 | off | 0.3 | +51 (score +20 and +51) |
| H1 | `rad_vet/time2_gtime_slope=true` | 0.3 | +51 |
| H0h | off | 0.15 | +102 (only if H1 - H0 in L_top is 0.9 to 2 %, to see which side of the dt limit H1 is) |
Stop H1 early (STOP file) if dt < 45 s for 5 cycles or the first stage fallback: that is already a FAIL.

## Return
- per cycle (hst + log): time, dt, L_top, L_mid, L_int (L_bot), CW, stage fallback / BE / dt_retry counts, s/cyc; the `SDO` lines
  (K1K2<0 share, max|ef/es-1|, cells > 0.1 / > 1, ef<=0) of both arms; GPU count, md5, .prov.
- L_top / L_mid / L_int of H1 against H0 at matched TIME (interpolated), at +5 +10 +20 +51, and max over the run; also vs the old-prod
  hst O30 of the Raven record (Orion stag2/hgtime/hst/O30.hst, scorer stag2/hgtime/lt.py) for orientation. Orion scores the hst.

## PASS / FAIL, fixed before the run
Record: release f, FHD2C: dt 90 throughout, L_top max 2.56 % vs O30 (IA9 2.69 %); allowed 2.92 % (old's Raven p99); f = 1/g FHD2: dt 27-51.
**FAIL (any one):**
1. dt: H1's median dt over cycles 15-51 < 60 s, or < 0.8 x H0's, or any cycle < 45 s while H0's is not.
2. Robustness: H1 has more CLEAN WARNING, stage fallbacks, BE redos or dt_retry than H0; any FATAL; any `ef<=0` cell in H1.
3. L_top: |H1 - H0| > 2 % at +20 or at +51 (matched time), or H1 outside old O30 by more than 2.92 % where H0 is inside.
4. Cost: s per simulated second > 1.05 x H0's.
**PASS** = none of these. Between 0.9 and 2 % in L_top is "run H0h", and H1 must then lie on the same side as H0h (toward the dt limit).
Not a criterion: the K1K2<0 share and max|ef/es-1| (no calibrated value; they explain a FAIL).
Hint for a FAIL: FHD2's drop came at ~14 cycles = ~1260 s, close to the 1170 s rotation block; report whether H1's drop aligns with a block
boundary (a coincidence in timing, not a shown cause).

-- orion (coordinator)
