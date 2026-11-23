# NOTE 2026-11-23i orion -> Caltech: USER decision: adopt IA11; restart bsg / m35mms / t9l62 FRESH from t = 0 on IA11

USER decision (11-23 ~09:50 PST): "adopt IA11 and restart caltech runs on it".
Restart all three FRESH from t = 0 (original ICs, same as IA9): no restart from the IA9 rst.

## 1. Stop the IA9 runs
- Stop bsg / m35mms / t9l62 (IA9) cleanly. Use a STOP file if your launcher has one. Otherwise wait until the run's latest periodic rst is complete,
  then cancel only its step. AthenaK has no signal-triggered final dump. Leave holder 7333072 alone.
- Keep their hst + logs as the IA9 reference: copy them into caltech_results/PROD_IA9_1123/<run>/ with a final PROD table row.
- Delete their rst except the latest one each. Post the kept rst names, t and md5.

## 2. Binary: GitHub branch `caltech-ia11-bd9553537` (head bd9553537)
- Build it exactly like IA9: CUDA 12.9 + g++ 13, Kokkos_ARCH_HOPPER90=On, -ffp-contract=off. Do not use CUDA 12.2 (the ray-closure miscompile).
- It is caltech-ia9-77551e823 plus 6 commits (src/rad_m1 + the he_star_m1 hst):
  6f99d454e clean_w_fse_tau, 5222d861d hst L_midE/L_intE, 957bbf7cb fwb F_ray of each solve's E, b9cbf7e3a fse_tau default on,
  08113c834 fwb_stage key-off fix, bd9553537 time2_gtime_mom.
- Post the build md5.

## 3. Tiny BSG gate (smokechk/runcpu.sh, close/bsgs.athinput, 1 GPU, as before)
Our IA11 binary on 1 A100 (bd9553537, CALQ recipe: mesh/meshblock nx2 = nx3 = 32, nlim 2):

    pass 0 9.534e-02 | pass 1 3.686e-11 | pass 2 1.026e-14 -> converged after 3 passes
    E 4.1527620104e+45 -> 4.1602693712e+45, F/c -> 2.9006e+42 1.1804e+39 -2.0678e+33

- E and F/c are the same as IA9 to every printed digit. The solution is the same and converges faster: pass 1 3.7e-11 instead of IA9's 1.885e-06.
- The faster pass 1 comes only from the fwb default flip. The same IA11 binary with clean_stag_fwb_stage = clean_stag_fwb_t2 = true in the
  gate deck gives IA9's numbers exactly (1.885e-06, 4 passes). bsgs.athinput names neither key, so IA11's default false applies.
- New reference for gate2.sh: pass 1 = 3.686e-11, F/c = 2.9006e42. Run it before the launch and post it.

## 4. Decks: `docs/handover/caltech_ia11/IA11_{bsg,m35mms,t9l62}.athinput` (md5 in MD5SUMS)
- Each is your IA9 deck (caltech_ia9/IA9_*.athinput) with only these changes, all in <rad_m1>:
  - clean_stag_fwb_stage / clean_stag_fwb_t2 = false. IA9 had true; IA11's code default is false, but the decks name them.
  - clean_w_fse_tau = true and time2_gtime_mom = true added (both are the IA11 defaults, named in the deck).
  - Plus a 2-line comment header.
- Everything else is identical: original ICs, cfl 0.3 in <time>, ray rotation off, GL 7x18, all keys in the deck, NOTHING on the command line.
- As for IA9, replace only the file paths (he_ic_file, he_opac_table, he_planck_table, + he_ic_linhao_file for bsg) with your copies of the
  same files. Change nothing else.
- Diff against IA9 (bsg; m35mms / t9l62 are the same, at lines 204-260 / 207-263):

```
0a1,2
> # IA11 2026-11-23 (USER: adopt IA11): GitHub caltech-ia11-bd9553537 = clean-stag2-interim-c bd9553537. = IA9_bsg.athinput (77551e823) with ONLY:
> #   rad_m1/clean_stag_fwb_stage = false, clean_stag_fwb_t2 = false (IA10 FDIP: thin rows F of each solve's own E); + clean_w_fse_tau = true (EPSEFF); + time2_gtime_mom = true (CALQ sdirk2 wall momentum)
240a243
> time2_gtime_mom = true   # IA11 (CALQ: sdirk2 stage-2 flux momentum f = 1/g)
243a247
> clean_w_fse_tau = true   # IA11 (EPSEFF: fse eps over the cell optical depth)
291,292c295,296
< clean_stag_fwb_stage             = true         # UNISET IA9 (all three stars)
< clean_stag_fwb_t2                = true         # UNISET IA9 (all three stars)
---
> clean_stag_fwb_stage             = false        # IA11 (FDIP; IA9 had true)
> clean_stag_fwb_t2                = false        # IA11 (FDIP; IA9 had true)
```

## 5. Launch
Use the same layout, holder and launcher as your 23d: 1 node = 4 H200 per run, 4 ranks, holder 7333072 (ends 11-26 05:21).
Use launch_ia9.sh with the IA11 binary and decks. New run dirs are fine (e.g. prod_ia11_1123/<run>), and so are the same rst/bin cadences.

## 6. Expected first cycles (our A100 smoke of exactly these decks)
Raven A100, bd9553537 (md5 35b66b73), 200 cycles, fresh. The only change is the meshblock split in x2/x3 for 40 GB GPUs: bsg / m35mms 8 blocks,
t9l62 32 blocks on 24 GPUs. Compared with the IA9 and IA10 smokes of the same layout (IA10 = IA11 without time2_gtime_mom):

| star | arm | rows | CW | FATAL | fallbacks | dt_retry | MCAP | s/cycle (10->end) | dt first / min / last | L_top/L_in end | L_bot/L_in end | never-read of the 4 keys |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bsg | IA9 | 201 | 0 | 0 | 0.0 | 0 | 1 | 0.932 | 165.46 / 164.68 / 168.83 | 1.0547 | 1.5731 | none |
| bsg | IA10 | 201 | 0 | 0 | 0.0 | 0 | 2 | 0.913 | 165.46 / 164.68 / 168.83 | 1.0547 | 1.5731 | none |
| bsg | IA11 | 201 | 0 | 0 | 0.0 | 0 | 2 | 0.909 | 165.46 / 164.64 / 168.72 | 1.0547 | 1.5341 | none |
| m35mms | IA9 | 201 | 0 | 0 | 0.0 | 0 | 0 | 0.972 | 12.008 / 12.008 / 12.008 | 1.0011 | 0.9993 | none |
| m35mms | IA10 | 201 | 0 | 0 | 0.0 | 0 | 0 | 0.948 | 12.008 / 12.008 / 12.008 | 1.0011 | 0.9993 | none |
| m35mms | IA11 | 201 | 0 | 0 | 0.0 | 0 | 0 | 0.945 | 12.008 / 12.008 / 12.008 | 1.0011 | 0.9993 | none |
| t9l62 | IA9 | 201 | 0 | 0 | 0.0 | 0 | 253 | 1.453 | 96.852 / 96.845 / 96.845 | 1.0022 | 0.8650 | none |
| t9l62 | IA10 | 201 | 0 | 0 | 0.0 | 0 | 228 | 1.393 | 96.852 / 96.845 / 96.845 | 1.0022 | 0.8650 | none |
| t9l62 | IA11 | 201 | 0 | 0 | 0.0 | 0 | 228 | 1.393 | 96.852 / 96.851 / 96.851 | 1.0022 | 0.9387 | none |

- In all three: 0 CW, 0 FATAL, 0 stage fallbacks, 0 dt_retry. Cycle-0 dt = IA9's (bsg 165.46, m35mms 12.008, t9l62 96.85).
- 14 keys are never read, the same set as your IA9 runs. None of the 4 changed keys is in that list. All three decks are sdirk2, so
  time2_gtime_mom is read; the time2_tableau line prints no ", time2_gtime_mom = false".
- bsg: L_bot/L_in at cycle 200 is 1.534 (IA9 1.573). This is the gtime_mom wall-momentum change during the start-up transient (see 7c).
  m35mms is unchanged. t9l62: L_bot/L_in at cycle 200 is 0.939, against 0.865 for IA9 and IA10. This is the wall fix; IA9 ran low.
  t9l62 cycle-0 dt is 96.852 and stays at 96.851 (IA9 96.845). s/cycle and MCAP count are the same as IA10.
- Your H200 numbers should match these to round-off in dt and in L_top/L_in, for the shared rows.

## 7. What changed vs IA9, and why
a. **fwb thin flux from each solve's own E** (clean_stag_fwb_stage = false, fwb_t2 = false; 957bbf7cb + 08113c834).
   The thin rows' F is the rays' flux factor times the E they leave: c (h+ E_L + h- E_R) of E^{n+1}, the old production's form.
   IA9 used the stage-old E. Its stage slopes put an O(h K) face-momentum dipole at HG's opaque photosphere, a flux dipole at fixed radius.
   This change removes it. HG L_top 2.66 % (F ratio 1.002 vs old 1.007); RSG 0.355 %; R136 0.102 %.
b. **fse optical-depth eps** (clean_w_fse_tau = true; 6f99d454e + b9cbf7e3a).
   fse's scattering isotropises only over the cell's optical depth: eps_eff = 1 - (1 - eps)(1 - exp(-tau_d)), tau_d = chi V / A_f.
   A thin or empty cell keeps the C1 face. This clears the dark-thin-domain solver cap (ws66 test). Star windows = IA9.
c. **sdirk2 wall momentum fix** (time2_gtime_mom = true; bd9553537).
   Under sdirk2, the hydro's stage-2 flux state took the stage-1 exchange's momentum with the energy's stiffness factor. In stiff cells that left
   the flux velocity a (1 - g) dt share of the residual radiation force short. The result was an O(dt) inner-wall drift: IA9's t9l62 L_bot/L_in
   was 0.72 at cfl 0.3, against 0.88 dt-converged, and the wall v^2 was damped about 100x.
   Now momentum uses f = 1/g in every cell, with the matching kinetic-energy change in IEN.
   Full-res t9l62, t = 0 -> 3e5 s: L_bot/L_in on J3/J4 within 0.002 in every window. Your IA9 run is 0.03-0.08 low. Wall v is not damped.
   Energy closes to 6e-9; restart + dt_retry bitwise; HG / RSG / R136 neutral.
- New hst columns (5222d861d): user.hst gains L_midE and L_intE at the end (2 more columns). These are the energy flux through the L_mid / L_int
  faces. All earlier columns keep their positions.

## 8. Please post
Build md5, the tiny gate, the kept IA9 rst, then command lines, nodes and first cycles (as in 23d). Then run the same ~20-min PROD table poster for the IA11
runs (new NOTE, columns as in 23g: CW, fallbacks, DT_RETRY, MCAP, FATAL, dt, s/cycle, L_top/L_in, L_bot/L_in, the last 3 hst rows), with the
hst copied to caltech_results/PROD_IA11_1123/<run>/. If anything FATALs or warns at startup, stop and post before retrying.

-- orion (coordinator)
