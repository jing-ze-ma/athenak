# TASK raven -> Caltech: He giant LARGE domain (x1max 700 Rsun), FRESH start from the MESA IC, FIN4, 4 H200

From the Raven coordinator (worker HGL). USER 10-14: fresh start from the MESA IC on a larger radial domain (no remap);
it runs on Caltech (4 H200). Raven hegR / Viper hegV He giant chains are cancelled. Reply to
NOTE-2026-10-14-caltech-hegiant-rmax.md. **Your current 200-Rsun chain (hgF4w + 4513229/30) and hgGuard keep running
until the guard stops them; do not cancel them for this TASK.**

## 0. Design (details: Raven /raven/ptmp/jinma/hegiant_large_1014/NOTE.md)
- **x1max = 700 Rsun** (4.8699e13 cm). Your numbers: R_ph 74.6/82.4/88.2/88.5 Rsun at 10.75/12.75/15.75/16.25 d,
  r(<rho>=1e-12) 98 Rsun at 16.25 d, ~1.2 Rsun/d and slowing. Conservative (no slowing) 1.2 Rsun/d gives r12 ~ 520 Rsun and
  R_ph ~ 510 Rsun at 1 yr; a log-decelerating fit gives ~200 Rsun. The rmax guard at 0.75 x1max = **525 Rsun** therefore covers
  ~1 yr in the conservative case. 500 Rsun would trip the guard (375) at ~250 d; 700 costs only +4 % cells vs 500
  (the cells sit in the fine band, not in the coarse top).
- Caveat (energy): the base heating (6.74e6 Lsun = L_MESA(3 Rsun)) equals the envelope's total energy above 3 Rsun
  (gas + rad - grav = -1.7e47 erg) after ~77 d. If the envelope then unbinds, the surface can accelerate (tens of km/s) and
  no affordable domain holds it for a year; the guard stops the run when r12 or R_ph reaches 525 Rsun.
- **Grid nx1 = 1697** (was 897), same StretchRPoly family and the same target law as fit2_897 (fit_giant3.py): dr/r inside
  100 Rsun = the current grid within 0.8-1.3x (60-95 Rsun fine band 580 vs 589 cells); the band law dr = H_atm (r/61.78)^2/3
  continues from 95 to 525 Rsun (dr/r 0.0010 at 100, 0.0016 at 150, 0.0020 at 200, 0.0029 at 300, 0.0037 at 400,
  0.0071 at 525; old grid 0.0031 at 100, 0.015 at 150, 0.019 at 200), then coarser to 0.0099 at 700. Monotone,
  max neighbour ratio 1.096. Code faces = fit faces to 7e-11 (he_grid_dump).
- **Sponge r0 = 630 Rsun** (0.9 x1max, E1 rule); rate 1e-3, he_sponge_dmax 1e-17 (> floor 1e-20: the floor gas is damped
  everywhere), outflow top (no inflow), marshak flux top: unchanged.
- **tfloor 1.2e3 K, eos_logt_min 3.0** (E1 rule "tfloor below T_rad at x1max": the code's floor gas sits at 1.4e3 K at
  690 Rsun after ic_rad; 3.0e3 K would clamp it beyond ~235 Rsun). Table nodes stay aligned with the old 3.4.
- **IC** ic_giant_own_x800.txt (md5 ac69934e) = your ic_giant_own.txt (ddc72692) BYTE-IDENTICAL to 230 Rsun + floor rows
  230-800 Rsun (rho 1e-20, T 4000 K, F_r = L_ph/4 pi r^2) from the same builder (make_ic_giant_mesa.py --rtop 800 Rsun).
- Everything else = the FIN4 deck hegiant-1014/hegiant_fin4.athinput.in (1b4d0592). tlim 1e9 (long-tlim rule).

## 1. Build
- Same binary as your current He giant: fin4-1014 @ 5c504e49 (or your d505d21f build, md5 1440c053: bitwise to 5c504e49).

## 2. Files (docs/handover/hegiant-large-1014/; md5)
- `hegiant_large.athinput.in` **005be1fc** -> `sed "s#@FILES@#<files dir>#g"`; change no key.
- IC: `cat ic_x800_head.txt <files>/ic_giant_own.txt ic_x800_tail.txt > <files>/ic_giant_own_x800.txt`;
  check md5 **ac69934e7d17951d708c4e917f95a46a** (head 077ddfd6, tail 242826ad, ddc72692 as before). Tables unchanged
  (rosseland c1ef063c, planck b699c60d).
- Reference: fit_giant3.py bdab6e24, fit3_1697.out d96515f6, faces3_1697.txt ce5b43a6, drr_table.txt 52fa9490,
  bg_large.py eada8978 (= hegiant-1014/bg.py with bands split at 200 and 373 Rsun).

## 3. Layout
- 1 node x 4 H200, 4 ranks, **1 MeshBlock 1697 x 64 x 64 per GPU** (mesh 1697 x 128 x 128), your CUDA_VISIBLE_DEVICES
  wrapper, --cpu-bind=cores, **export UCX_TLS=^xpmem**. No lowmem key (Viper MI300A: 82 GB per 1697x64x64 block;
  your 897 run had 85.7 GB for 2 blocks of 897 -> expect ~80-85 GB of 141 GB per H200).

## 4. Gate (one job; report each item)
A. Smoke C: fresh, `time/nlim=24 time/ndiag=1 output1/dt=1e-6 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output6/dt=1e30
   output5/dt=240`. Report rc, FATAL/NaN/DIVERGED/NON-CONVERGED, Picard mean/max, dt, s/cycle (3-24), GPU peak memory, the
   startup lines `tfloor_kelvin = 1200 K` and `sponge rate = 0.001 above r = 4.38291e+13`, and the vet_gd halo line
   ("band-clamped lateral reads ... 0").
B. Restart R: new srun `-r <C>/rst/hegiant.00001.rst time/nlim=24` -> hst rows after cycle 12 + final dumps BITWISE vs C.
   Not bitwise -> STOP and report.
C. Background BG (~2 h wall): fresh, deck + `problem/he_seed=0.0 output2/dt=1000 output3/dt=1e30 output6/dt=1e30 output5/dt=1e30`,
   `-t 01:55:00`. Then `BASE=hegiant LOGF='run*.log' VISPY=<repo>/vis/python python3 docs/handover/hegiant-large-1014/bg_large.py
   4.30368e12 4.38291e13 <BG dir> 1`. Pass if: bands 0-1.5 R_ph match the E0 reference at equal t (hegiant-1014/bg_E0_ref_raven.txt;
   E0 at 1.0e4 s: 0.99-1.01 +2.69 km/s / 0.510, 1.03-1.10 -0.17 / 0.035, 1.1-1.5 -0.15 / 0); the floor bands 1.5 R_ph-r0 and
   >r0 have |<v_r>| <= 0.2 km/s and max|drho/rho| <= 0.01 (sponge-terminal infall, g/rate: 0.03 km/s at 100 Rsun -> 0.001 at
   690 in the smoke); Mdot_top ~ 0; NC 0; budget L_in + P_esrc = 6.755e6 Lsun (x 12 wedge).
   (Viper bg 12206134 on the half-phi wedge is queued; Raven will send its numbers if it runs first.)

## 5. Production chain (start after A+B pass and BG passes; no GO file needed unless BG fails)
- Run dir e.g. `/resnick/groups/carnegie_poc/jingze/hegiant_large_1014/run`; fresh t = 0; your link pattern (run.cfg read at
  job start; FRESH if no rst else newest valid rst + rst_info.py last_time keys; STOP on rc != 0 / FATAL / NaN / DIVERGED /
  NON-CONVERGED > 0; prune keep 3).
- **tlim 1e9** in the deck and run.cfg TLIM=1e9; **7-d links** (QoS normal), athena -t = 7 d minus your margin; link 1 +
  **2 afterany backups**.
- **rmax guard**: your hegiant_1014/rmax_guard.sh with **RLIM=525** (0.75 x1max), pointed at the new run dir (STOP + scancel of
  the new chain's jobs only).
- Cost: s/cycle ~ nx1 -> ~0.9 s/cycle on 4 H200 (your 897: 0.48); dt 20 s (unchanged: set by the inner lateral cells; smoke
  dt 20.0) -> ~6.7e5 cycles = ~1.3e7 s (~155 d) of star time per 7-d link; 1 yr ~ 2.4 links (~17 d wall).
- Storage: dumps 0.25 d (hydro_w ~0.56 GB, m1 ~0.45, m1_vet ~0.45 per dump) ~1.45 GB per 0.25 d -> ~2.1 TB per star-year;
  rst ~19 GB each (keep 3). Keep the deck cadences unless storage forbids; tell Raven if you change any.
- Page: add the large run to the He giant page (https://claude.ai/artifact/9g8xcCyBCjdzwkBbALaiqt) as a new section
  ("x1max 700 Rsun, fresh"), R_ph(t) and r12(t) on the same axes as the 200-Rsun run.

## 6. Report: NOTE-2026-10-14-caltech-hegiant-large.md
- build md5, deck/IC md5 after substitution; A/B/C results; chain job ids; first-link s/cycle, dt, memory.
