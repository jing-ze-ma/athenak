# TASK raven -> Caltech: E1W = 2x-LATERAL narrow wedge started by REMAPPING the running e1s state; 1 H200, 7-d links

From the Raven coordinator (worker RMP, e1remap_1014). Files: `docs/handover/e1w-1014/`. Report NOTE:
`docs/handover/NOTE-2026-10-14-caltech-e1w.md`.

USER (10-14): run the 2x-lateral narrow-wedge AG Car E1 (pi/6 x pi/6, theta 5pi/12..7pi/12, phi 0..pi/6, 536 x 128 x 128 = ONE
MeshBlock on 1 H200) starting from a REMAP of the running e1s state, not from t = 0. Source = the NEWEST e1s restart file at the time
you run this TASK (whatever t; post-ramp is fine). Everything else (deck physics, FIN4 + dt_retry 5b7e1bc5 nofma numerics, tlim 1e9)
is as e1s.

## Rules
- Do NOT touch e1s (agcE1s 4513330-32 and its prod/ dir: read only; you only COPY one rst out of it), BSG, He giant. Cancel nothing.
- New dir `/resnick/groups/carnegie_poc/jingze/e1w_1014/{src_rst,remap,scripts,gate,prod}`.
- Never put an ABSENT key on the command line (every key you override below exists in the generated deck).

## 1. Build (patch on 5b7e1bc5, nofma, as your e1s binary)
- `e1w-1014/e1remap_b06ccb37.patch` (md5 05a25cec09ddf97683490add855c96a1) = `git diff 5b7e1bc5 b06ccb37` (branch e1remap-1014;
  ONLY src/pgen/he_star_m1.cpp). In a COPY of your 5b7e1bc5 source tree: `patch -p1 < e1remap_b06ccb37.patch`;
  then `md5sum src/pgen/he_star_m1.cpp` must be **59e3334b331dc45f4e1e01cff5b642d7** (5b7e1bc5's is 65eeb9ba...).
- Build exactly as `athena_gpu_he_star_m1_5b7e1bc5_nofma` (same modules, cmake flags incl. nofma) ->
  `builds/bin/athena_gpu_he_star_m1_b06ccb37_nofma`; record its md5 (differs from 6252e58e; that is expected).
- What the patch does: new key `problem/he_ic_remap_file` (read only when NAMED; absent = the old code path, bitwise). When named on a
  FRESH start: reads the remap file (header checked against <mesh>; FATAL unless `<time>/start_time` equals the parent's t bitwise),
  takes the parent's pgen restart block (he_ic_equil a_ref -> Phi_eff exactly as a restart; Picard window stripped), skips ic_rad /
  he_ic_equil / seed, loads u0, M1 u0, f0x1 (ghosts incl.). The rest is a fresh start: hesdirk2 first step BE (time2_start = be),
  implicit F0 faces re-derived from the remapped cell F, E, v (implicit_f0_init = auto), empty rotavg/csk/pred history, dt_retry fresh;
  w_mlt(t) and the vet_gd rotation (physical clock) continue from t. Restarts of the new run are ordinary restarts.

## 2. Remap (python3 + numpy; ~2-5 min, < 16 GB RAM; on a CPU or the gate node)
1. `R=$(ls -t <e1s>/prod/rst/*.rst <e1s>/keep_rst/*.rst | head -1)` (newest by time); `cp -p $R e1w_1014/src_rst/`; md5sum both
   (must match; if e1s was writing it, take the previous one). Record name, t (`remap.txt` prints it), md5.
2. Template deck `e1w-1014/agcA_E1W_fin4dtr_s3_tlim1e9.athinput` md5 **09403077791531cde03fd146e6145ee9** (the narrow-wedge deck:
   680dadb0 + x2min/x2max/x3max + meshblock 536x128x128 + seed k band; Caltech file paths).
3. `python3 e1w-1014/remap_e1w.py e1w_1014/src_rst/<rst> agcA_E1W_fin4dtr_s3_tlim1e9.athinput e1w_1014/remap --name e1w > remap/remap.txt`
   (script md5 **5c10eeb42d42827326465d22f81480f6**; bin_convert.py beside it not needed here). It writes `remap/e1w.remap`
   (778,720,240 bytes) + `remap/e1w.athinput` (= template + `<time> start_time = t`,
   `<problem> he_ic_remap_file = <abs path>`, and every `<outputN>` `last_time` / `file_number` aligned with e1s's schedule and
   numbering). `chmod a-w` both. Paste remap.txt into the NOTE. Expected (section 4): Phi spread <= ~1e-12, per-shell
   conservation of rho/ET/E <= ~1e-13 (eint is NOT conserved by design: the KE of the sub-cell velocity, capped at 1 % of the
   parent eint, comes out of it), momentum per parent cell ~1e-14, blend conservation <= ~1e-12, seam metric before ~2-15 ->
   after 0, MeshBlock-border check ~0.5-1 (4 parent blocks; >> 1 would mean misplaced blocks), no `SystemExit`. Runtime ~2 min.
   - Method: central pi/6 x pi/6 sub-wedge of e1s (64 x 64 parent cells; theta centred on pi/2; phi offset chosen automatically as the
     least seam mismatch, printed), radial grid copied, conservative 2x lateral prolongation (volume centroids, MC + Barth-Jespersen
     limited) of rho, eint (rst slab), E, F, f0x1; velocity = parent mean + limited deviation (momentum exact per parent cell),
     its sub-cell KE capped at 1 % of eint and taken from eint; ET rebuilt with the parent's rho*Phi (ET exact per parent cell);
     periodic seams by a conservative pairwise blend over W = 8 child cells each side (q_L' = (1-a) q_L + a q_R, a = 0.5 (1 - d/8)):
     seam jump 0, every shell sum exact.

## 3. Gate (one job, 1 x 1 H200, -t 2:00; `export UCX_TLS=^xpmem`)
- Launch `e1w_1014/launch1.sh`: `exec srun --mpi=pmix -n 1 -c 8 --cpu-bind=cores /resnick/groups/carnegie_poc/jingze/e1_1013/e1wrap.sh "$@"`.
- G0 key-absent bitwise: the TEMPLATE deck (no remap key), fresh, `time/nlim=4 time/ndiag=1 output1/dt=1 output2/dt=1 output3/dt=1
  output4/dt=1e30 output5/dt=1e30`, once with 6252e58e and once with the new binary: every hst row + bin BITWISE and the final rst
  (`rst/agcar3d.00001.rst`) byte-identical (`cmp`). (Init ~1-2 min each.)
- G1 smoke + restart: `smoke_e1.sh` (copy e1s_1014/scripts as for e1s) with launch1.sh and `remap/e1w.athinput`: C (60 from the
  remap) vs L (30 + 30 hook restart) [+ Coff]; `smoke_gate.sh`: CHAINCMP / CMPLAST / every hst row BITWISE -> VERIFIED. Check the init
  log lines `he_ic_remap_file ... t = <t>`, `he_ic_equil: a_ref ... (restart trailer)`, `he_ic_remap_file: u0, M1 u0, f0x1 loaded`,
  `w(t = <t>) = <w>`, `band 102 cells (needed 102)`, end `band-clamped ... 0`; dt (expect ~401-410 s); peak GB.
- G2 transient: fresh from `remap/e1w.athinput` (directly with launch1.sh, in gate/G2), `time/nlim=300 time/ndiag=10
  output2/dt=4.1e3 output2/last_time=-1 output3/dt=1e30` (~1.2e5 s, ~6-7 min, 30 hydro_w dumps ~5 GB). Run
  `python3 e1w-1014/ana_remap.py gate/G2 8` (bin_convert.py beside it; agcar3d.x1grid.txt is written by the run): per dump rms v_r in
  the SEAM band (< 3 cells from a seam) and the INTERIOR (>= 12 cells), at 0.4-0.95 and 1.0-1.2 R_ph, max|v_r|, KE_lat, KE_rad,
  rho seam metric. Also user.hst row 0 vs the e1s user.hst row at the same t: L_bot, L_top x 3.8638 within ~1 % of e1s, w_mlt equal (the in-loop hst rows of e1s report w of the step just taken, i.e. one dt earlier: ~5e-5 offset is expected).
  PASS (as in the A536 test, section 4): no jump at start (KE_lat, KE_rad change < 5 % over the first 10 cycles), seam-band v_r
  <= ~1.2x interior and not growing faster than interior, rho seam metric < 1.5, FATAL/NaN/NON-CONVERGED 0, DT_RETRY 0 (a recovered
  retry: report, not a STOP). Paste the ana_remap table (every 5th row) into the NOTE.
- Write `prod/GO_e1` ("USER 10-14 e1w remap, TASK-2026-10-14-caltech-e1w, gate PASS <date>") only if G0, G1, G2 pass. Else STOP, report.

## 4. Raven/Viper validation (b06ccb37; logs in `e1w-1014/test_logs/`)
| test | where | result |
|---|---|---|
| key ABSENT, e1w deck fresh 4 cyc, 1 MI300A: 053bd996 (= 5b7e1bc5 numerics) vs b06ccb37 | Viper 12210570 | hst + 12 bins BITWISE, final rst BYTE-IDENTICAL |
| key ABSENT, 268x8x8 fresh 4 cyc, CPU: 5b7e1bc5 vs b06ccb37 | Raven 31172615 | hst + 12 bins BITWISE, final rst BYTE-IDENTICAL |
| FULL SIZE: e1s deck fresh 6 cyc on 2 MI300A (4 MeshBlocks 536x64x64, rst 5.1 GB, t 4364 s) -> remap -> e1w 536x128x128 | Viper 12210849 (final script; 12210569 same with the earlier script) | remap: MB-border check 0.84/0.69 (rho), Phi spread 4.5e-13, shell conservation rho 1.0e-14, ET 1.2e-14, E 1.6e-14, momentum 2.8e-14, blend 9.5e-15 (m1 5e-13); C 20 cyc on 1 MI300A: rc 0, FATAL/NaN/DT_RETRY 0, a_ref from the trailer, band 102 (needed 102), band-clamped 0, dt 401.25, Picard 3.0, init 27 s, ~1.9-2.0 s/cyc, peak 70.5 GB; L 10 + restart 10: CMPRUNS BITWISE (hst 21/21, 12 bins) |
| CONVECTIVE state: A536 (536x32x32 pi/3, seed 1e-2, rst t 9.0005e6 = infall onset, vrms 8-9 km/s) -> 536x32x32 pi/6 | Raven CPU 31172851, Viper 12210823 | remap: shell conservation rho/ET 2.7e-13, E 2.2e-14, momentum 7e-15, blend <= 1.6e-12; 225 / 138,752 parent cells KE-capped; seam metric [phi, theta] rho 2.6/4.0, v_r 3.5/4.3 -> 0/0. C 30 vs L 15 + restart 15 (CPU): CMPRUNS BITWISE (hst 18/18, 10 bins). 400 cycles (3.8e5 s) on 1 MI300A: rc 0, 0 DT_RETRY, Picard 7.5 (parent 9.5) |

A536 transient (child G vs the parent continued from the same rst, P; rms v_r km/s, KE erg; child KE x 3.8638 for the parent wedge):
| t (s) | G env seam / int | G surf seam / int | G KE_lat, KE_rad | G rho seam [phi, th] | P env | P KE_lat, KE_rad |
|---|---|---|---|---|---|---|
| 9.0005e6 (remap) | 4.51 / 9.00 | 8.68 / 9.09 | 1.74e43, 2.45e43 | 0, 0 | 7.49 | 1.05e44, 1.37e44 |
| 9.008e6 (8 cyc) | 4.51 / 9.00 | 8.72 / 9.13 | 1.75e43, 2.45e43 | 0.13, 0.11 | 7.49 | 1.05e44, 1.37e44 |
| 9.10e6 (110 cyc) | 4.46 / 9.12 | 9.34 / 9.74 | 1.82e43, 2.50e43 | 0.28, 0.27 | 7.52 | 1.06e44, 1.37e44 |
| 9.376e6 (400 cyc) | 4.75 / 9.46 | 11.5 / 12.2 | 2.15e43, 2.87e43 | 0.44, 0.30 | 7.53 (9.15e6) | 1.07e44, 1.37e44 |
- No start-up jump (KE change 0.4 % in the first 8 cycles); the seam band starts at ~0.5x the interior v_r in the envelope (the blend
  averages two uncorrelated sides) and stays there; rho seam jumps stay below the typical adjacent jump (metric < 0.5).
- Child KE_lat grows ~23 % over 3.8e5 s while the parent's is flat: the doubled lateral resolution resolving more lateral motion
  (smooth growth, no jump). The central sub-wedge holds ~2/3 of the per-sr KE of the full parent wedge at this t.
- user.hst row 0: L_bot x 3.8638 = 2.8516e38 vs parent 2.8520e38; L_top x 3.8638 = 2.762e38 vs 2.751e38 (local sub-wedge, 0.4 %).
- Reduced-radius children (nx1 268) cannot be tested: the vet_gd band is geometric (radial cell size) and needs 38 > 32 / 73 > 64
  lateral cells (FATAL clamped reads at cycle 0). The production geometry 536 x 128 x 128 needs 102 <= 128, independent of the state.

## 5. Production (`e1w_1014/prod`)
- run.cfg: BIN = the b06ccb37 nofma binary (its md5), IN = `remap/e1w.athinput` (the generated deck; it is a FRESH start at t),
  KEYS empty, TLIM=1e9, STOP_NC=1, E1C_GATE=1, E1C_VERIFIED=gate/smoke/VERIFIED, LAUNCH=e1w_1014/launch1.sh, E1C_DIR=e1w_1014/scripts,
  E1C_MARGIN=1200, `export UCX_TLS=^xpmem`. prod/rst empty at the first link (link 1 starts from the remap deck; later links restart).
- Links `agcE1w`: `-A carnegie_poc -p gpu --qos=normal --exclude=hpc-sm-01-09 --ntasks=1 -c 8 --gres=gpu:nvidia_h200:1 --mem=120G
  -t 7-00:00:00` (QoS long = 1 job/user, held by e1s), first link + 2 backups `afterany`, same walltime.
- Keep `remap/` and `src_rst/` until the first prod rst exists (then src_rst may go; keep remap.txt).
- Push a STARTED note: jobs, node, t_start, first s/cycle (expect ~1.2 on H200), GB, dt, band/clamp lines.

## 6. Reports (side by side with e1s at matching times)
At t_start + 5e5, + 2e6, then every ~5e6 s: one table per time with an e1w row and the e1s row at the same t
(both from `pulse_ana.py run` + user.hst): vrV, vrM, vr105, KEe, rsurf, L_top/L0, KE_int (e1w x 3.8638 = Omega_e1s/Omega_e1w,
1.04720/0.27103 sr; give raw too), w_mlt, Picard mean, DT_RETRY, s/cycle, dt. Bins share e1s's numbering (dump n at ~n x 1e5 s),
so compare same-numbered dumps. Shell means (vr*, d*, rsurf) and ratios compare directly; wedge totals need x 3.8638.
Also the e1w-1014 cells-per-H_p / dominant lateral wavelength / v_rms / L_kin comparison at 0.80, 0.86, 0.92 R_ph (IC expectation
cells per H_p e1s 8.9/6.7/3.3, e1w 17.8/13.4/6.6). Push prod_pulse.txt + agcar3d.user.hst to `docs/handover/e1w-1014/caltech_out/`.
Flag at once: infall (vrV or vrM < -5 km/s), rsurf > 1.15, DT_RETRY bursts, Picard mean > 20, band-clamped reads > 0.

## NOTE skeleton
build md5 (patch, he_star_m1.cpp, binary) | src rst (name, t, md5) + remap.txt | G0/G1 BITWISE lines | G2 transient table |
STARTED | evolution tables (e1w vs e1s).
