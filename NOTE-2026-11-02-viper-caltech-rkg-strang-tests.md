# NOTE 2026-11-02 viper -> Caltech (USER): please run the RKG Lie-vs-Strang GPU tests (branch rkgstrang-test-1102 @ e204131b2); viper GPUs are crowded

User: "move tests to Caltech, it has free GPUs". Thanks for b9fad64ad / a20499535 (fin51 = aeafae97 noted).

## What it is
AthenaK runs the RKG STS of the resistive EMF (Ohmic + ambipolar) ONCE over the full dt AFTER the RK stages
(driver.cpp): a Lie split, first order in dt, although RKG2 is second order. Athena++ (RKL2, its default) is Strang:
0.5*dt before the main integrator and 0.5*dt after (main.cpp sts_loc op_split_before/after; sts_task_list.cpp:751-753).
New key <mhd>/rkg_strang (default false = unchanged): dt/2 before the RK stages, dt/2 after; each half takes its
stage count from dt/2; eta refreshed before the second half; rkg_bflx = mean of the halves; FATAL with <rad_m1>.
Commits on aeafae97: ee383da2c (option), e204131b2 (output variable mhd_eta = eta_O + eta_A, diagnostic).
CPU convergence test (advected resistive wave: Lie / Strang / no-STS orders) runs on viper.

## Please run (H200, your r1m production keys incl. use_rkg_sts, rkg_lean, rkg_dt_safety 0.5; any DEVELOPED MHD
state, e.g. the off_0p50.rst you used for mhdspeed; same binary e204131b2 for all arms)
Gates:
- G1: rkg_strang absent, 40 cycles: e204131b2 vs aeafae97 -> rst data after <par_end> + dhj.mhd.hst IDENTICAL.
- G2: rkg_strang = true: 20 + 20 vs 40 cycles -> rst data IDENTICAL.
Accuracy arms (each ~0.02 rotation from the same restart; rkg_strang must be in the deck, a command-line key that is
not in the input FATALs):
- lie (current), strang, lie4 (time/cfl_number / 4 = reference), strang4 (cfl / 4), seed (lie + problem/seed = 1,
  seed_amp = 1e-14, seed_restart = true; = the noise floor), strangB (strang repeated: GPU noise check).
Per pressure band (p < 1e-6, 1e-6..1e-3, 1e-3..1, > 1 bar), volume-weighted:
- RMS over cells of (arm - lie4) for B_r, B_zonal, B_merid, v (3 comps), T; and relative difference of the band
  totals of E_B and of the resistive heating q;
- the same for seed - lie (noise), so each arm is quoted as a multiple of the noise.
Verdict wanted: is |lie - lie4| above the noise where it counts (p > 1e-6 bar), and does strang cut it (expect
error ratio lie/lie4 ~ 4x per cfl/4 if first order, strang ~16x)? Also ms/cycle lie vs strang, interleaved x2 on a
clean GPU (expect strang ~1.4x the RKG stages).

## Report
NOTE on bsg-files-1009 with the gate verdicts, the band table and timing. Nothing to merge: rkg_strang stays
default off until the user decides.
