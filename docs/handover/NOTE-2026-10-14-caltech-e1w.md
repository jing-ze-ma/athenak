# NOTE Caltech -> Raven: E1W gate (TASK-2026-10-14-caltech-e1w, 4d8efb28) -- G0 / G1 / G2 PASS, GO_e1 written, links NOT submitted

Caltech test-runner, 2026-10-14. Dir `/resnick/groups/carnegie_poc/jingze/e1w_1014/{src_rst,remap,scripts,gate,prod}`.
e1s (agcE1s 4513330-32) untouched: one rst copied out, nothing else read-write. Production links are held for the user's go
(the coordinator submits `prod/submit_links.sh`).

## 1. Build
- Source: copy of the persistent 5b7e1bc5 tree `builds/inc_gpu_he_star_m1_nofma/src` (he_star_m1.cpp 65eeb9bafb0ff300f68c4e3a627f6b21)
  -> `builds/inc_gpu_he_star_m1_e1w_b06ccb37/src`; `patch -p1 < e1remap_b06ccb37.patch` (md5 05a25cec09ddf97683490add855c96a1);
  he_star_m1.cpp -> **59e3334b331dc45f4e1e01cff5b642d7** (as required).
- Fresh CMake configure with the 5b7e1bc5_nofma flags (gcc 13.2 + cuda 12.9 + hpcx 2.17.1, nvcc_wrapper, HOPPER90, MPI, Release,
  CMAKE_CXX_FLAGS=-ffp-contract=off, Kokkos d8e9af03; compile line and Kokkos options identical to the 5b7e1bc5 build).
  expansion job 4527789 (16 min).
- Binary `builds/bin/athena_gpu_he_star_m1_b06ccb37_nofma` md5 **6d4215c032c6620f4e13c94d48e03af6** (5b7e1bc5: 6252e58e).

## 2. Remap
- Source rst: `e1s_1014/prod/rst/agcar3d.00045.rst` (newest; written 13:39, copied 13:55, size and mtime unchanged across the copy), t = **22500363.63636364**
  (cycle 31564), md5 **4f38c89a0c69606dc5ae50eb4e0d0fe7** (original and copy identical). 4,506,865,759 bytes.
- Template deck md5 09403077791531cde03fd146e6145ee9, remap_e1w.py md5 5c10eeb42d42827326465d22f81480f6. expansion job 4527825,
  37 s, MaxRSS 2.6 GB. `remap/e1w.remap` 778,720,240 bytes md5 26e23d31a4770813a0ca329e43e5faca; `remap/e1w.athinput` md5
  ee3977115602c0f62ef4a388cddee815 (start_time 22500363.63636364; hst last_time 2.249e7 / file 2250; bins + log 2.24e7 / 225;
  rst 2.2e7 / 45). Both chmod a-w. No SystemExit.
- Versus section 4 expectations: Phi spread 4.5e-13 (ok), shell conservation rho 5.0e-14 / ET 5.2e-14 / E 5.7e-15 (ok),
  momentum per parent cell 6.0e-14 / 3.3e-14 / 4.3e-14 (Viper 2.8e-14), blend <= 8.9e-14 (ok), seam metric 5.4-7.9 -> 0 (ok),
  KE cap active in 1561 / 2,220,032 parent cells. MeshBlock-border check theta 0.97 / phi 1.29 (rho), E phi 1.43: slightly above
  the 0.5-1 range, not >> 1 (this state is convective, parent seam metric itself ~5.7).

remap.txt:
```
parent rst /resnick/groups/carnegie_poc/jingze/e1w_1014/src_rst/agcar3d.00045.rst size 4506865759
  t = 22500363.63636364  dt = 727.27272727247328  ncycle 31564  nmb 4  mesh (536, 128, 128)  meshblock (536, 64, 64)  ng 3
  blocks ['PGENST01', 'M1TIME2A', 'EINTRST1', 'M1RSSIG1', 'M1CSK1S1', 'M1ROTAV1']  nhydro 5  f0x2/f0x3 (1, 1)  counts {'npred': 0, 'nt2': 13, 'neint': 1, 'nctr': 0, 'nrss': 5, 'nvgd': 0, 'ncsk': 3, 'nra': 16}
  parent theta [1.0471975511966, 2.0943951023932] 128, phi [0, 1.0471975511966] 128; child theta [1.30899693899575, 1.83259571459405] 128, phi [0, 0.523598775598299] 128
  sub-wedge: parent theta cells 32..95 (64), phi 64 cells
  MeshBlock-border check rho  rms jump across borders / elsewhere: theta 0.970 phi 1.288
  MeshBlock-border check eint rms jump across borders / elsewhere: theta 0.960 phi 1.279
  MeshBlock-border check E    rms jump across borders / elsewhere: theta 0.914 phi 1.430
  rho*Phi from the parent (ET - KE - eint): Phi lateral spread (max rel, active) 4.535e-13; min eint 1.5319e-05, min rho 2.7409e-17
  phi0 auto: parent phi cell offset 42 (seam mismatch / mean adjacent jump 3.488; k0=0: 4.062; median 3.862)
  velocity slopes: KE added by the sub-cell velocity structure taken from eint, capped at 0.01 eint_p: max 1.000e-02 of eint_p; alpha < 1 in 1561 of 2220032 parent cells (min alpha 0.284)
  momentum m1 per parent cell: max |sum children - parent| / max|m| 5.999e-14
  momentum m2 per parent cell: max |sum children - parent| / max|m| 3.270e-14
  momentum m3 per parent cell: max |sum children - parent| / max|m| 4.325e-14
  prolongation conservation per shell (max rel over all radial cells) rho   5.041e-14
  prolongation conservation per shell (max rel over all radial cells) ET    5.163e-14
  prolongation conservation per shell (max rel over all radial cells) E     5.693e-15
  prolongation conservation per shell (max rel over all radial cells) eint  1.185e-03
  prolongation conservation per shell (max rel over all radial cells) f0x1  5.735e-15
  blend conservation per shell rho   max rel 7.643e-15
  blend conservation per shell ET    max rel 8.055e-15
  blend conservation per shell E     max rel 9.360e-15
  blend conservation per shell m1    max rel 8.885e-14
  blend conservation per shell f0x1  max rel 1.003e-14
  after blend: min eint 1.5541e-05 (min eint/parent-min 1.015), min rho 6.6465e-17, min E 1.0000e-04
  seam metric = rms(seam jump)/rms(adjacent jump), active radial cells: [phi seam, theta seam] before -> after (W=8)
    rho  [5.659, 5.664] -> [0.00e+00, 0.00e+00]
    ET   [5.377, 5.202] -> [0.00e+00, 0.00e+00]
    E    [5.803, 6.056] -> [0.00e+00, 0.00e+00]
    m1   [6.220, 6.272] -> [0.00e+00, 0.00e+00]
    v_r  [7.946, 6.199] -> [0.00e+00, 0.00e+00]
    rho phi-adjacent jumps: rms in blend band 4.353e-10, rms interior 5.039e-10
    v_r phi-adjacent jumps: rms in blend band 6.409e+04, rms interior 6.681e+04
  child wedge totals (lateral-weighted sums, radial factor omitted): rho 1.2771077311e-04 ET 4.4020622208e+10 E 3.4810536779e+10
  wrote /resnick/groups/carnegie_poc/jingze/e1w_1014/remap/e1w.remap (778720240 bytes) md5 26e23d31a4770813a0ca329e43e5faca; pgen block 4328 bytes
  wrote deck /resnick/groups/carnegie_poc/jingze/e1w_1014/remap/e1w.athinput md5 ee3977115602c0f62ef4a388cddee815
```

## 3. Gate (GPU job 4527945 on hpc-sm-02-10, 1 x 1 H200, `launch1.sh` = srun -n 1 -c 8 e1wrap.sh, UCX_TLS=^xpmem; G1 rerun 4532039)

### G0 key absent, template deck fresh 4 cycles: 5b7e1bc5 (6252e58e) vs b06ccb37 (6d4215c0) -- PASS
```
hst agcar3d.hydro.hst: rows A 5 B 5 common 5 differing 0
hst agcar3d.user.hst: rows A 5 B 5 common 5 differing 0
bins compared 12 differing 0
CMPRUNS BITWISE
cmp: all 12 bin files IDENTICAL; rst/agcar3d.00000.rst IDENTICAL; rst/agcar3d.00001.rst IDENTICAL (byte-identical)
```

### G1 smoke + restart from remap/e1w.athinput -- PASS (on the 2nd run; 1st run failed on a harness cadence artefact)
- 1st run (job 4527945, kept in `gate/smoke_dt4e3`): smoke_e1.sh test cadence hst 4e3 / bins 1.5e4 is NOT on the remap deck's
  output grid (last_time 2.249e7 / 2.24e7 are not multiples of 4e3 / 1.5e4); the hook re-aligns the restarted link to
  floor(t/dt)*dt, so C and L wrote rows / dumps at different times: chaincmp "missing 3, extra 2, differing 0" -> GATE FAIL.
  All common rows were identical and CMPLAST BITWISE: a schedule mismatch, not a data difference. Production cadence (hst 1e4,
  bins 1e5, rst 5e5) IS on the deck grid, so the production hook alignment is consistent.
- 2nd run (job 4532039): `scripts/smoke_e1w.sh` = smoke_e1.sh with the one test-key line changed to output1/2/3 dt = 5e3
  (5e3 divides both last_time values), nothing else changed. smoke_gate.sh:
```
hst agcar3d.hydro.hst: C rows 6, chain rows 6, missing 0 [], extra 0 [], differing 0
hst agcar3d.user.hst: C rows 6, chain rows 6, missing 0 [], extra 0 [], differing 0
bin: C dumps 50, chain dumps 50, extra in chain 0 []
CHAINCMP BITWISE
CMPLAST BITWISE
HST agcar3d.hydro.hst BITWISE
HST agcar3d.user.hst BITWISE
GATE PASS 2026-10-14_15:19:26    (Coff: CMPLAST BITWISE, IPC-OFF HST hydro/user BITWISE)
```
- Init lines (C, Coff; L link 1 the same):
  `he_ic_remap_file .../e1w.remap: t = 2.25004e+07, pgen block 4328 bytes (remapped fresh start)`;
  `he_ic_equil: a_ref (cells is .. ie) and Phi_eff from the fresh start's equilibrium (restart trailer)`;
  `he_ic_remap_file: u0, M1 u0, f0x1 loaded (1 MeshBlock(s) on rank 0); no ic_rad, no he_ic_equil, no seed`;
  `w(t = 2.25004e+07) = 0.0172725` (L link 2: w(t = 2.25124e+07) = 0.0170288);
  `vet_gd: exact per-shell lateral halo, band 102 cells (needed 102)`; end `band-clamped lateral reads ... 0`.
  FATAL 0, NaN 0, NON-CONVERGED (final counter) 0, Picard mean 27.96 (C/Coff 60 cyc), 21.0 (L link 2).
- **DT_RETRY: one recovered retry, deterministic, in every run (C, L link 2, Coff, G2) at cycle 40:**
  `Picard NON-CONVERGED after 200 passes: resid=4.17e-07 (tol 1e-08)` -> `DT_RETRY cycle=40 t=2.251651e+07 dt 410.26 -> 208.86
  reason=NON-CONVERGED cut 1` -> `DT_RETRY recovered cycle=54 ... cap off`. Reproduced bitwise across the restart. Reported, not a STOP.
- dt 401.25 (cycle 10) -> 410.26 / 411.19 (expected ~401-410).

### G2 transient, fresh from remap/e1w.athinput, 300 cycles (t 2.25004e7 -> 2.26207e7, 1.2e5 s) -- PASS
- rc 0, FATAL 0, NaN 0, NON-CONVERGED (final) 0, band 102 (needed 102), band-clamped 0, DT_RETRY 1 recovered (the cycle-40 one
  above), Picard mean 20.8 (max 152), hesdirk2 BE steps 46, stage fallbacks 45.
- `agcar3d.x1grid.txt` is not written by the G2 run (the smoke C run did write one); the identical-grid file from gate/smoke/C was
  copied into gate/G2 for ana_remap.py (same mesh x1 lines). ana_remap.py run on expansion (job 4532018).
- ana_remap table (W = 8; every 5th dump + last; rms v_r km/s, KE erg, raw e1w wedge):
```
#  t(s)        cyc  vr_env_seam vr_env_int vr_sur_seam vr_sur_int  max|vr|  KE_lat     KE_rad     rho_seam_phi rho_seam_th
2.250040e+07      0     4.2329     6.9255     2.8168     2.7993    38.654 2.027e+43 2.897e+43    0.000    0.000
2.252100e+07     57     4.2126     6.9432     2.9994     2.8679    39.117 2.027e+43 2.914e+43    0.335    0.367
2.254150e+07    107     4.2095     6.9644     3.1748     2.9929    39.787 2.022e+43 2.934e+43    0.278    0.270
2.256210e+07    157     4.2350     6.9826     3.3272     3.1256    40.662 2.023e+43 2.954e+43    0.280    0.283
2.258260e+07    207     4.2672     6.9994     3.4466     3.2516    41.325 2.028e+43 2.974e+43    0.296    0.347
2.260310e+07    257     4.2933     7.0160     3.5335     3.3813    41.913 2.034e+43 2.995e+43    0.307    0.349
2.262070e+07    300     4.3120     7.0304     3.5788     3.4870    42.441 2.042e+43 3.011e+43    0.327    0.346
```
  (2nd dump, 11 cycles: KE_lat 2.030e43, KE_rad 2.900e43.)
- Criteria: KE_lat +0.15 %, KE_rad +0.10 % over the first 11 cycles (< 5 %); envelope seam band 0.61x interior and flat (as A536);
  surface seam/interior 1.006 -> 1.026 (<= 1.2x), seam +27 % vs interior +25 % over 300 cycles (same rate; both rising smoothly);
  rho seam metric max 0.57 (< 1.5), settles 0.27-0.35.
- user.hst row 0 (t = 22500363.6363636) vs e1s user.hst row at the same t (x 3.8638):
  L_bot 7.38531e37 x 3.8638 = 2.85354e38 vs e1s 2.85398e38 (-0.015 %); L_top 3.27073e37 x 3.8638 = 1.26374e38 vs 1.25424e38
  (+0.76 %); w_mlt 0.01727254 vs 0.01728731 (e1s row = w of the step just taken; offset -1.5e-5); KE_int x 3.8638 = 0.924 of e1s.

### Cost
- s/cycle on 1 H200 (G2, ndiag 10, a hydro_w dump every ~10 cycles): 2.06-2.7, median ~2.3 (cycles 60-300); smoke C (ndiag 1,
  dumps every ~12 cycles): (239.1-81.8)/50 = 3.1. Higher than the ~1.2 the TASK expects; production cadence (dumps every ~240
  cycles) should come in lower; the STARTED note will give the production number.
- Init ~60-80 s. Peak GPU memory 65.7 GB (nvidia-smi, 2-s sampling; G0+G2 job 64.8 GB). dt ~401-411 s.
- Note: a fresh start from the remap deck writes an init rst `rst/agcar3d.00045.rst` (1.32 GB, smaller than the 4.2 GB link-end
  rst; G0 shows the same for a t = 0 start, 3.1 vs 4.75 GB).

## GO_e1 and production (prepared, NOT submitted)
- `prod/GO_e1`: "USER 10-14 e1w remap, TASK-2026-10-14-caltech-e1w, gate PASS 2026-10-14". VERIFIED = gate/smoke/VERIFIED.
- `prod/run.cfg`: BIN b06ccb37_nofma, BIN_MD5 6d4215c032c6620f4e13c94d48e03af6, IN remap/e1w.athinput, KEYS "", TLIM 1e9, STOP_NC 1,
  E1C_GATE 1, E1C_VERIFIED gate/smoke/VERIFIED, LAUNCH e1w_1014/launch1.sh, E1C_DIR e1w_1014/scripts, E1C_MARGIN 1200,
  export UCX_TLS=^xpmem. prod/rst empty.
- `prod/submit_links.sh`: agcE1w, -A carnegie_poc -p gpu --qos=normal --exclude=hpc-sm-01-09 --ntasks=1 -c 8
  --gres=gpu:nvidia_h200:1 --mem=120G -t 7-00:00:00, link 1 + 2 backups afterany, link script e1w_1014/scripts/link_e1.sh.
- Keep remap/ and src_rst/ until the first prod rst exists.
