# RESULTS 2026-09-30: DeltaAI WASP-121b 1x rot 300 -> resistive MHD, max_eta scan (bbot 3 G / 10 G)

All 13 DeltaAI arms finished and analysed. Task: rt-integration `docs/handover/TASK-2026-09-29-deltaai-w121-mhd-etamax.md`.
Only the viper arms (b10_e14_f13, b10_e13_f16, b3_e13_f16; held on viper) are still `TBD(<arm>)`; fill them with
`collect_arm.sh <arm> [OHMIC_DT]` (this dir) and `docs/handover/scripts/w121/mhd_eta/ana_eta.py runs/<arm>`.
Per-arm `ana_eta.txt` in `ana_eta/`; job scripts and history (`run_arm_deltaai.sub`, `rotate_interactive.sh`,
`rotation.log`, `done.txt`) in this dir. ana_eta.py uses the DeltaAI dME/dt window fix (+-0.25 rot, hst cadence 0.1 rot).

## Setup
- Run root `/work/nvme/bivj/jma20/w121_mhd_0929` (`runs/<arm>/`: `run.<jobid>.log`, `dhj.mhd.hst`, `dhj.user.hst`,
  `ana_eta.txt`; job history `rotation.log` / `rotation.txt`; pre-fix crashed outputs in `failed_mpireq/`).
- Binary `/u/jma20/ATHENAK/builds/bin/athena_dhj_dev_2dfe4e95ff1e`, md5 **7ba630ea1f10349619b448a2c71b1f73**
  (source rt-integration 793e03c3 = 6e7eeb3e + seam-flux request fix 2dfe4e95). Every finished arm's job .out
  records this md5.
- 1 node GH200, ghx4-interactive, one arm at a time. **2 ranks = 2 GH200 per arm** (12 MeshBlocks each), not the
  4 of `run_arm.sub`. The noise twin b3_e12_f13_r2 ran on **1 rank = 1 GH200** (log: `Number of parallel ranks = 1`;
  rotation.txt NRANK 1).
- Package as shipped: lhlld, `ohmic_resistivity = eos`, `cs_lowbeta_fallback 0.5`, `<mhd>/dfloor 1e-13` (user 09-29),
  `hlld_bx_zero_tol` unset (1e-4), cfl 0.3, fresh `-i` start from remap.dat (identity remap of w1x rst 00600),
  tlim = rot 302 on the command line (the input's `tlim` is rot 300; do not read it). P_rot 110153.5 s.
- Arm status:
  - finished + analysed (13): b3_e0/e11/e12/e12_r2/e13/e13_sts/e14_sts_f13, b10_e0/e11/e12/e13/e13_sts/e14_sts_f13.
    All rc 0, `Terminating on time limit` at rot 302; single link except the two below.
  - b10_e13_sts_f13: 2 links. Job 3276005 stopped on the WALL limit at t = 3.325791e7 = rot 301.923, **96.1 % of tlim**
    (not 99.6 % as noted in memory); -r finish job 3277988 (~9.9k cycles, `finish_b10_e13_sts.sh`) reached tlim.
  - b10_e14_sts_f13: 3 links. Job 3278345 stopped on the WALL limit at t = 3.325012e7 = rot 301.853 (92.6 %),
    -r job 3279924 on the WALL limit at t = 3.326306e7 = rot 301.970 (98.5 %), -r job 3280176 reached tlim.
  - viper (binary ae767d20): b10_e14_f13, b10_e13_f16, b3_e13_f16; HELD by the user 09-30 01:25 (SESSION-2026-09-30-viper.md).
  - **b3_e14_f13: not run anywhere.** Only `runs/b3_e14_f13/mhd.athinput` exists here; it was never in
    rotation.txt / queue.txt. viper's job 12028522 was cancelled by the user 09-29 23:40 (split NOTE: "b3_e14_sts_f13
    covers it"). So the 3 G e14 explicit cost is estimated below, not measured.
- Failed pre-fix links (binary 952abbfe, no seam fix): b3_e0/e11/e12 aborted in MPI_Isend (Cray MPICH "Internal
  error") at cycle 8100 = the seam-flux request leak (NOTE-2026-09-29-deltaai-seam-mpi-leak.md); the old
  b3_e12_f13_r2 (1 rank) was cancelled by hand, not crashed. All four reran from scratch on 2dfe4e95.

## Caveat: all arms predate the lhllc/lhlld radial-face fix
- rt-integration aa9d51b7 (NOTE-2026-09-30-lhllc-x1-phi.md) adds `<mhd>/lhlld_x1_phi_min`, default **1 for dhj**:
  full HLLD velocity-jump term on x1 (radial) faces. Our binary (793e03c3, an ancestor of aa9d51b7) has phi ~
  chi(2-chi) ~ 3e-3 on radial faces in the deep envelope, i.e. a **2-cell radial velocity checkerboard is
  undamped** (hydro 1x: odd-even share .9-.96 at 0.1-1e4 bar, ~97 % of the radial KE). The rot-300 source state
  (hydro lhllc) carries it too.
- What it means here:
  - All arms share the same solver and start, so the **relative** cap comparisons (dt, cost, trends with cap) stand.
  - **Absolute** numbers that are sensitive to grid-scale radial motion are suspect: E_mag, Q_ohm, the
    grid-scale current and eta_num, and deep (>10 bar) flow (viper: deep v_h +6 % with the fix, still growing).
    A grid-scale radial v shears B at the grid scale, so the 1-10 bar field / heating is the part most at risk.
  - Jet / day-night / dt: viper saw <1 % jet change and unchanged dt with the fix (hydro; MHD smoke key 1:
    dt unchanged, radial KE 10x lower after 100 cycles).
  - Any restart of these rst files with a post-aa9d51b7 binary silently switches to key 1 (dhj default applies on
    restart). Continuing an arm with the new binary is therefore a solver change mid-run.
  - Before production: rerun the chosen-cap arm and its e0 reference per bbot with the fixed binary (~15-25 min
    each on 2 GH200) and check that E_mag, Q_ohm/L_abs and the cap-sweep fractions move by less than the cap effect.

## Acceptance criterion (TASK, verbatim in substance)
The acceptable cap is the LARGEST cap that meets all of:
- (a) stable: no FATAL, NaN or dt collapse; div B at round-off.
- (b) wall/rot <= ~2x the eta ~ 0 arm (same bbot, dfloor), or affordable with STS (compare STS twins).
- (c) energy-containing scales cap-independent vs the next-higher cap and the e0 arm at equal time: jet max,
  u_eq at 1e-3/0.1/1 bar, day-night contrast at 0.1 bar, Q_ohm/L_abs, at-cap t_drag vs t_adv; the cap-sweep
  spurious-drag fraction should be ~0 for p > 1e-4 bar. Only p > 1e-6 bar counts. Differences below the noise
  twin spread are noise.

## Per-arm cost and stability (run logs; collect_arm.sh)
Stability = counts of NOT-CONVERGED / FATAL / nan / "dt collapse" over all run.*.log of the arm. dt from the
100-cycle `elapsed=` lines (cycle-weighted); "late" = last 0.5 rot (rot 301.5-302, after the winding transient).
Wall/rot includes startup (6-10 s). Limiter: explicit Ohmic dt = cfl/6 min(dx^2)/eta_cap = 17.5 s (1e12),
1.78 -> 1.755 s (1e13), 0.18 s (1e14) on these states; "Ohmic" = dt within 3 % of it (flat dt), else CFL.

| arm | GPUs | NC/FATAL/nan/dtcoll | dt mean / min [s] | dt late [s] | limiter | wall min/rot (2 rot / late) | ms/cycle | wall/e0 (2 rot / late) |
|---|---|---|---|---|---|---|---|---|
| b3_e0_f13      | 2 | 0/0/0/0 | 4.116 / 2.913 | 3.47 | CFL (9.64 s at rot 300 -> 3.1-3.5 s after winding) | 6.5 / 7.6 | 14.46 | 1 / 1 |
| b3_e11_f13     | 2 | 0/0/0/0 | 3.414 / 1.960 | 2.25 | CFL (Ohmic 175 s) | 7.7 / 11.6 | 14.32 | 1.18 / 1.53 |
| b3_e12_f13     | 2 | 0/0/0/0 | 3.562 / 2.402 | 3.12 | CFL (Ohmic 17.5 s) | 7.4 / 8.4 | 14.30 | 1.14 / 1.11 |
| b3_e12_f13_r2  | **1** | 0/0/0/0 | 3.569 / 2.256 | 3.02 | CFL | 9.8 / 11.5 | 18.97 | (1 GPU) |
| b3_e13_f13     | 2 | 0/0/0/0 | 1.760 / 1.755 | 1.755 | **Ohmic 100 %** of samples | 14.7 / 14.8 | 14.13 | **2.26 / 1.95** |
| b3_e13_sts_f13 | 2 | 0/0/0/0 | 2.742 / 2.066 | 2.56 | CFL; RKG s=3 in 97 % (s=4: 2.4 %, 5: 0.7 %) | 17.1 / 18.3 | 25.55 | 2.63 / 2.41 |
| b3_e14_f13     | - | not run | est. 0.18 (Ohmic) | - | Ohmic | est. ~140 (late) | ~14 | est. ~19x |
| b3_e14_sts_f13 | 2 | 0/0/0/0 | 2.370 / 1.873 | 2.12 | CFL (Ohmic 0.18 s); RKG s=7 88 % (8: 3.9 %, 9: 2.1 %, 10: 4.4 %, up to 16) | 32.7 / 35.4 | 42.16 | **5.03 / 4.66** |
| b10_e0_f13     | 2 | 0/0/0/0 | 1.731 / 1.055 | 1.27 | CFL (3.09 s at rot 300) | 14.9 / 20.4 | 14.09 | 1 / 1 |
| b10_e11_f13    | 2 | 0/0/0/0 | 1.551 / 0.936 | 1.08 | CFL | 16.8 / 24.0 | 14.20 | 1.13 / 1.18 |
| b10_e12_f13    | 2 | 0/0/0/0 | 1.384 / 0.911 | 0.97 | CFL | 18.7 / 26.9 | 14.13 | 1.26 / 1.32 |
| b10_e13_f13    | 2 | 0/0/0/0 | 1.095 / 0.848 | 0.95 | Ohmic only in the first ~4 % of cycles, then CFL | 23.8 / 27.4 | 14.18 | **1.60 / 1.34** |
| b10_e13_sts_f13 | 2 | 0/0/0/0 (2 links: 3276005 wall-stop at rot 301.92 + 3277988 -r finish) | 1.113 / 0.827 | 0.94 | CFL; RKG s=3 99.6 % | 42.0 / 49.8 (sum of both links, excl. restart startup) | 25.44 | 2.82 / 2.44 |
| b10_e14_f13 (viper) | - | TBD(b10_e14_f13) | est. 0.18 (Ohmic) | - | Ohmic | est. ~7x e0 | TBD | TBD(b10_e14_f13) |
| b10_e14_sts_f13 | 2 | 0/0/0/0 (3 links: 3278345 + 3279924 wall-stops, 3280176 -r finish) | 0.971 / 0.626 | 0.76 | CFL (Ohmic 0.18 s); RKG s=5 71 %, s=4 25 % (s=4 only after rot 301.85), s=6-11 4 % | 60.9 / 72.7 (sum of 3 links, excl. restart startup) | 32.21 (link 1 32.6, links 2-3 28.8) | **4.09 / 3.56** |
| b10_e13_f16 / b3_e13_f16 (viper) | - | TBD(b10_e13_f16) / TBD(b3_e13_f16) | v_A-limited 0.14 / 0.46 s est. | | | | | |

Notes on cost:
- ms/cycle is cap-independent without STS (14.1-14.5 on 2 GH200); cost differences are all dt.
- STS costs 25.5 ms/cycle (1.8x) with s=3 almost always (sts-lean-0930 cuts this ~10 %, not used here). At 1e13
  STS gains dt only at 3 G (2.56 vs 1.755 s late) and **nothing at 10 G** (CFL 0.95 s < Ohmic 1.755 s), so
  **STS is slower than explicit at 1e13 for both bbot** (3 G: 17.1 vs 14.7 min/rot; 10 G: 41.6 vs 23.8).
- 1 GH200 vs 2 GH200 (twin): 18.97 vs 14.30 ms/cycle, i.e. 2 GPUs give only 1.33x at 12 MeshBlocks/GPU.
- b3_e14_f13 estimate: late Ohmic dt 0.18 s vs e0 3.47 s at the same ms/cycle -> ~19x e0. b10_e14_f13 estimate:
  0.18 vs 1.27 s -> ~7x. Both fail (b) by far; e14 only via STS.

## Per-bbot comparison at equal time, rot 302 (ana_eta.txt; rot 300 = common start)
Only the first and last bins (rot 300, 302) were analysed (ana_eta.py default); **rot 301 is not in these
files** (needs `ANA_BINS=all`, written into a copy of the arm dir since ana_eta.py writes ARMDIR/ana_eta.txt).
rot 300 (all arms): jet max 14207 m/s at 4.2e-4 bar; u_eq 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar)
day/night 2039/1481, contrast 558 K.

### 3 G (bbot 0.846 code), dfloor 1e-13
| arm | jet max [m/s] (p bar) | u_eq 1e-3/0.1/1 bar [m/s] | contrast 0.1 bar [K] | E_mag [erg] | \|B\|max | Q_ohm/L_abs | eta_num(est) | max\|divB\|dr/\|B\| |
|---|---|---|---|---|---|---|---|---|
| e0        | 13073 (4.2e-4) | 8576/5071/2351 | 566 | 2.356e33 | 363 | 2.4e-13 | -2.8e10 | 6.3e-12 |
| e11       | 12767 (7.5e-4) | 8757/5046/2353 | 569 | 1.718e33 | 194 | 3.04e-3 | -2.9e11 | 7.0e-12 |
| e12       | 12881 (1.3e-3) | 9216/5066/2354 | 564 | 1.028e33 | 118 | 3.29e-3 | -7.3e11 | 7.6e-12 |
| e12 r2 (1 GPU) | 12763 (1.3e-3) | 9135/5066/2354 | 564 | 1.023e33 | 118 | 3.30e-3 | -7.4e11 | 4.9e-12 |
| e13       | 13606 (7.5e-4) | 9576/5089/2362 | 564 | 3.808e32 | 101 | 2.54e-3 | -1.0e12 | 1.0e-11 |
| e13 sts   | 13692 (7.5e-4) | 9585/5088/2362 | 564 | 3.801e32 | 101 | 2.48e-3 | -1.0e12 | 2.7e-11 |
| e14 sts   | 13819 (7.5e-4) | 9584/5086/2363 | 564 | 2.482e32 | 95.8 | 1.88e-3 | -9.4e11 | 1.6e-10 |

### 10 G (bbot 2.821 code), dfloor 1e-13
| arm | jet max [m/s] (p bar) | u_eq 1e-3/0.1/1 bar [m/s] | contrast 0.1 bar [K] | E_mag [erg] | \|B\|max | Q_ohm/L_abs | eta_num(est) | max\|divB\|dr/\|B\| |
|---|---|---|---|---|---|---|---|---|
| e0      | 11395 (7.5e-4) | 8079/4923/2314 | 588 | 7.827e33 | 520 | 1.0e-12 | -2.0e10 | 9.3e-12 |
| e11     | 11398 (1.0e-3) | 8167/4946/2325 | 584 | 5.170e33 | 304 | 7.93e-3 | -2.9e11 | 9.5e-12 |
| e12     | 11197 (4.2e-3) | 8065/5016/2334 | 571 | 3.216e33 | 196 | 8.17e-3 | -6.5e11 | 6.6e-12 |
| e13     | 11676 (3.2e-3) | 8785/5092/2350 | 561 | 1.486e33 | 192 | 6.69e-3 | -7.5e11 | 1.5e-11 |
| e13 sts | 11607 (3.2e-3) | 8668/5094/2350 | 562 | 1.490e33 | 190 | 6.72e-3 | -7.5e11 | 1.1e-10 |
| e14 (viper arm, not run on DeltaAI) | TBD(b10_e14_f13) | TBD(b10_e14_f13) | TBD(b10_e14_f13) | TBD(b10_e14_f13) | TBD | TBD(b10_e14_f13) | TBD | TBD(b10_e14_f13) |
| e14 sts | 12785 (1.0e-3) | 9303/5113/2354 | 560 | 1.122e33 | 155 | 5.46e-3 | -7.2e11 | 2.0e-10 |

### Noise (twin b3_e12_f13 1 GPU vs 2 GPU, rot 302)
jet max 0.9 %, u_eq(1e-3) 0.9 %, u_eq(0.1, 1 bar) 0, contrast 0, E_mag 0.5 %, Q_ohm 0.3 %, dt mean 0.2 %,
cap-sweep and at-cap fractions <= 0.007. Two rotations are too short for the chaos to saturate: treat this as a
LOWER bound on the noise.

### f(eta at cap), at-cap median log10 t_drag/t_adv (capped / uncapped), f(t_drag < t_adv) at cap; rot 302
Bands 1e-6-1e-4 / 1e-4-1e-2 / 1e-2-3e-2 / 3e-2-0.1 / 0.1-0.3 / 0.3-1 / 1-10 bar ("f cap/unc f<").
| arm | 1e-6-1e-4 | 1e-4-1e-2 | 0.01-0.03 | 0.03-0.1 | 0.1-0.3 | 0.3-1 | 1-10 |
|---|---|---|---|---|---|---|---|
| b3 e11 | .10 -4.42/-1.14 1.00 | .19 -4.44/-1.21 1.00 | .48 -4.14/-1.48 1.00 | .57 -3.78/-1.53 1.00 | .82 -3.26/-1.88 1.00 | 1.00 -2.81/-1.64 1.00 | .99 -2.44/-1.74 .98 |
| b3 e12 | .07 -3.43/-0.55 1.00 | .16 -3.43/-0.95 1.00 | .36 -2.76/-0.50 1.00 | .40 -2.42/-0.48 1.00 | .43 -2.13/-0.60 .99 | .53 -1.80/-0.66 .96 | .26 -1.59/-1.05 .95 |
| b3 e13 | .07 -2.94/-2.27 .99 | .10 -1.94/-0.60 1.00 | .27 -1.09/+0.51 .97 | .26 -0.77/+1.00 .92 | .25 -0.46/+0.82 .83 | .20 -0.27/+0.48 .74 | .03 -0.50/-0.13 .77 |
| b3 e13 sts | .06 -2.98/-2.18 1.00 | .10 -1.92/-0.60 1.00 | .28 -1.07/+0.51 .97 | .26 -0.75/+1.00 .91 | .25 -0.44/+0.83 .83 | .20 -0.25/+0.49 .73 | .03 -0.47/-0.10 .76 |
| b3 e14 sts | .04 -1.57/-0.65 .32 | .06 -0.42/+0.53 .79 | .15 +0.36/+2.40 .17 | .17 +0.73/+2.29 .05 | .13 +0.93/+2.07 .04 | .06 +0.89/+1.43 .05 | .00 +0.58/+0.63 .04 |
| b10 e11 | .11 -4.58/-0.44 1.00 | .21 -4.63/-1.16 1.00 | .48 -4.52/-1.69 1.00 | .58 -4.34/-2.11 1.00 | .83 -4.04/-2.63 1.00 | 1.00 -3.64/-2.41 1.00 | .99 -3.32/-2.59 1.00 |
| b10 e12 | .09 -3.86/-1.22 1.00 | .17 -3.74/-1.07 1.00 | .36 -3.18/-1.07 1.00 | .40 -2.90/-1.11 1.00 | .44 -2.67/-1.32 1.00 | .55 -2.49/-1.56 1.00 | .27 -2.32/-1.80 1.00 |
| b10 e13 | .07 -3.15/-2.21 1.00 | .11 -2.40/-0.95 1.00 | .26 -1.54/-0.08 1.00 | .26 -1.26/+0.27 .99 | .25 -1.02/+0.04 .96 | .20 -0.88/-0.24 .93 | .03 -0.91/-0.62 .93 |
| b10 e13 sts | .07 -3.11/-2.16 1.00 | .11 -2.39/-0.96 1.00 | .26 -1.56/-0.06 1.00 | .26 -1.27/+0.27 .99 | .25 -1.02/+0.02 .96 | .20 -0.88/-0.24 .93 | .03 -0.91/-0.64 .93 |
| b10 e14 sts | .02 -2.23/-0.77 .70 | .06 -0.87/+0.18 .91 | .15 -0.17/+1.64 .56 | .16 +0.17/+1.83 .44 | .13 +0.23/+1.40 .36 | .06 +0.25/+0.85 .26 | .00 +0.10/+0.17 .40 |
| b10 e14 (viper arm, not run on DeltaAI) | TBD(b10_e14_f13) | | | | | | |
(e0 arms: f(cap) = 1 everywhere by construction; uncapped at-cap t_drag/t_adv = -5.5 .. -1.5 (3 G), -5.9 .. -2.5 (10 G).)

### Cap sweep on each arm's rot-302 state: spurious-drag volume fraction, p > 1e-4 bar
Bands 1e-4-1e-2 / 0.01-0.03 / 0.03-0.1 / 0.1-0.3 / 0.3-1 / 1-10 bar (1e-6-1e-4: 0.01-0.04 everywhere).
| state | cap 1e12 | cap 1e13 | cap 1e14 | cap 1e15 |
|---|---|---|---|---|
| b3 e0   | .054 .087 .087 .067 .035 .015 | .053 .086 .086 .064 .025 .002 | .052 .085 .084 .057 .016 0 | .050 .063 .048 .025 .001 0 |
| b3 e12  | .024 .121 .144 .138 .101 .022 | .024 .119 .139 .128 .067 .003 | .023 .112 .110 .061 .016 0 | .021 .046 .009 .004 0 0 |
| b3 e13  | .018 .159 .198 .202 .192 .036 | .018 .148 .170 .158 .104 .006 | .017 .083 .037 .015 .009 0 | .003 .002 0 0 0 0 |
| b3 e13 sts | .019 .159 .198 .203 .192 .036 | .019 .148 .170 .157 .103 .006 | .018 .078 .035 .014 .009 0 | .003 .001 0 0 0 0 |
| b3 e14 sts | .050 .220 .240 .234 .220 .049 | .048 .179 .140 .107 .060 .006 | .032 .023 .008 .005 .003 0 | .005 .001 0 0 0 0 |
| b10 e0  | .060 .094 .066 .014 .002 .002 | .060 .094 .066 .014 .002 0 | .060 .093 .065 .013 0 0 | .057 .086 .062 .012 0 0 |
| b10 e12 | .029 .091 .100 .070 .035 .004 | .029 .091 .100 .068 .031 .001 | .029 .090 .096 .056 .013 0 | .028 .067 .029 .003 0 0 |
| b10 e13 | .014 .114 .140 .126 .093 .013 | .013 .112 .135 .112 .062 .003 | .013 .102 .093 .038 .013 0 | .007 .006 .002 0 0 0 |
| b10 e13 sts | .014 .116 .140 .125 .092 .013 | .013 .112 .134 .112 .062 .003 | .013 .102 .097 .040 .013 0 | .007 .006 .002 0 0 0 |
| b10 e14 sts | .041 .168 .182 .170 .135 .021 | .040 .152 .143 .114 .070 .003 | .032 .063 .055 .036 .012 0 | .009 .005 .001 0 0 0 |
| b10 e14 (viper arm, not run on DeltaAI) | TBD(b10_e14_f13) | | | |
Ohmic dt of the sweep caps on these states: 17.46 / 1.755 / 0.181-0.185 / 0.019 s.

## Criteria
(a) Stable: **all 13 finished arms pass** (0 NOT-CONVERGED/FATAL/nan/dt collapse; rc 0) after 54k-227k cycles.
Max |div B| dr/|B|: explicit arms 4.9e-12 to 1.5e-11, no drift with cap; STS arms higher, 2.7e-11 (b3 e13 sts),
1.1e-10 (b10 e13 sts), 1.6e-10 (b3 e14 sts), 2.0e-10 (b10 e14 sts) - still tiny, but ~5-10x the explicit arm at the
same cap. The viper arms (b10_e14_f13, f16): TBD.

(b) Cost vs e0 (2-rot mean / late):
- 3 G: e11 1.18/1.53, e12 1.14/1.11, **e13 2.26/1.95 (explicit)**, e13 sts 2.63/2.41, e14 sts **5.03/4.66**
  (final; partial was ~4.6/4.7), e14 explicit est. ~19x. -> e13 explicit is at the ~2x limit (passes on the late rate, which
  is the one production pays); e14 fails (b3_e14_sts final 5.03x / 4.66x).
- 10 G: e11 1.13/1.18, e12 1.26/1.32, **e13 1.60/1.34 (explicit)**, e13 sts 2.82/2.44 (final, 2 links),
  e14 sts **4.09/3.56** (final, 3 links; est. was ~3x), e14 explicit TBD(b10_e14_f13) (viper arm, not run on DeltaAI) (est. ~7x).
  -> e13 explicit passes; e14 fails.
- At 10 G the extra cost of e11-e13 is NOT Ohmic (Ohmic dt 175/17.5 s, and 1.755 s > the late CFL dt): the
  CFL dt itself falls with the cap (1.27 -> 1.08 -> 0.97 -> 0.95 s late). See "Odd".

(c) Cap independence:
- Flow at p >= 0.1 bar and the day-night contrast are nearly cap-independent: 3 G e12 -> e13 u_eq(0.1 bar) +0.5 %,
  u_eq(1 bar) +0.3 %, contrast 0 K; e0 -> e13 within 0.5 % / 0.5 % / 2 K. 10 G e12 -> e13 +1.5 % / +0.7 % / -10 K;
  e0 -> e13 +3.4 % / +1.6 % / -27 K (-4.6 %). Twin noise there is ~0 (lower bound).
- The upper atmosphere is cap-dependent above noise: u_eq(1e-3 bar) rises monotonically with the cap
  (3 G: 8576/8757/9216/9576, +3.9 % e12 -> e13 vs noise 0.9 %; 10 G: 8079/8167/8065/8785, +8.9 %); jet max
  e12 -> e13 +5.6 % (3 G) and +4.3 % (10 G). At 1e14 (STS) the 3 G 1 mbar flow has saturated (u_eq(1e-3) 9584 vs 9585
  e13 sts, jet +0.9 %) while 10 G keeps rising (9303, +5.9 % vs e13 explicit, +7.3 % vs e13 sts; jet 12785, +9.5 % vs e13).
  10 G e13 -> e14 sts at >= 0.1 bar: +0.4 % / +0.2 % / -1 K. The field brakes the rot-300 jet (14207) by 4-10 % (3 G) and 18-21 % (10 G).
- E_mag is NOT converged in the cap: 3 G 2.36/1.72/1.03/0.38 e33 (e0/e11/e12/e13), 10 G 7.83/5.17/3.22/1.49 e33;
  each decade removes ~30-60 %, more at higher caps; e14 sts 0.25 e33 (3 G, -35 % vs e13 sts) and 1.12 e33 (10 G, -25 %). |B|max likewise. Q_ohm/L_abs is small at every cap (2.5-3.3e-3 at 3 G, 6.7-8.2e-3
  at 10 G) and flat e11-e13; it drops at e14 sts (1.88e-3 at 3 G, 5.46e-3 at 10 G).
- t_drag test FAILS at 1e13 for both bbot: 20-27 % of the 0.01-1 bar cells sit at the cap; in them the capped
  t_drag/t_adv is 10^-1.1 .. 10^-0.3 (3 G) / 10^-1.5 .. 10^-0.9 (10 G) while the uncapped value is 10^+0.5 .. 10^+1.0
  (3 G) / ~10^0 (10 G): the cap creates drag the physical eta would not.
  Cap sweep on the e13 states, 0.01-1 bar: 10-17 % spurious at 1e13, 1-10 % at 1e14, <= 0.6 % only at 1e15
  (target "~0 for p > 1e-4 bar"). The smoke-state expectation (0 at 1e13, "10 G moves one decade") does not hold
  for the wound field; both bbot need ~1e15 by this test, and on the e0/e12 states even 1e15 leaves up to 5-9 % at 1e-4-0.1 bar.
- t_drag test at 1e14 (e14 sts arms, own rot-302 state): 6-17 % of the 0.01-1 bar cells sit at the cap (both bbot); in
  them the capped t_drag/t_adv is 10^+0.4 .. 10^+0.9 (3 G) / 10^-0.2 .. 10^+0.3 (10 G), and f(t_drag < t_adv) at cap is
  4-17 % (3 G) / 26-56 % (10 G). Cap sweep at 1e14 on these states, 0.01-1 bar: 0.3-2.3 % spurious (3 G, below the
  1-10 % predicted from the e13 states), 1.2-6.3 % (10 G, inside it); 1e-4-1e-2 bar 3.2 % for both.
- Strictly, no affordable cap meets (c). 1e14 (measured, e14 sts arms) cuts the spurious fraction to 0.3-2.3 % (3 G) /
  1.2-6.3 % (10 G) at 0.01-1 bar, still not ~0.

## Provisional recommendation (DeltaAI arms all done; pending only the viper arms b10_e14_f13, f16)
- **3 G: max_eta = 1e13, explicit (use_rkg_sts = false). PROVISIONAL.** Largest cap passing (a) and (b)
  (1.95x late, 2.26x incl. winding). STS is 16 % slower at 1e13 (s=3 floor). (c) is met for winds >= 0.1 bar and
  the day-night contrast, not for E_mag, the 1 mbar flow, or the t_drag test. Upgrade to 1e14 + STS only if
  b3_e14_sts_f13 comes in at <= ~2x e0: final 5.03x (2 rot) / 4.66x (late) - not met.
- **10 G: max_eta = 1e13, explicit. PROVISIONAL.** (a) ok, (b) 1.34x late / 1.60x, and at late time the run is
  CFL-limited so the cap is free there; STS costs 1.8x more per rot at 1e13 (final: 42.0 vs 23.8 min/rot = 1.76x, late 49.8 vs 27.4 = 1.82x). 1e14 needs
  b10_e14_sts_f13 (cost; measured 4.09x / 3.56x late) and TBD(b10_e14_f13) (viper arm, not run on DeltaAI) (explicit reference for STS
  accuracy); est. was ~3x (STS) / ~7x (explicit) e0 -> fails (b) for STS (measured).
- Say explicitly to the user: with an affordable cap the magnetic energy (and anything that depends on the field
  amplitude, e.g. Ohmic heating per band, magnetic drag above 1e-2 bar) is set by the cap, not by the physical
  eta. The dynamical impact on the >= 0.1 bar winds is small (<~2 %) over 2 rotations; longer runs may differ.
- Recheck after the lhlld radial-face fix (caveat above) before production. The f16 check arms (viper, held)
  decide whether the dfloor 1e-13 choice biases the top: TBD(b3_e13_f16), TBD(b10_e13_f16).

## Odd / to check
- ana_eta.py "dt run" at the rot-302 bin is the hst dt of the LAST step, which run.sub clips to hit tlim
  (e.g. b3_e13: 0.875 s vs 1.755 s in the log), so the rot-302 "limiter" label is wrong for b3_e13 ("CFL/other";
  the log shows Ohmic 100 %). At rot 300 the STS arm is labelled CFL correctly only by accident. Use the log dt
  (table above), not the ana limiter line.
- CFL dt falls with the cap where Ohmic is not binding: 10 G late dt 1.27/1.08/0.97/0.95 s for e0/e11/e12/e13, and
  at 3 G e11 (late 2.25 s) < e12 (3.12) < e0 (3.47). A higher cap gives LESS E_mag, so this is not the bulk v_A.
  Unexplained; the CFL-limiting cell is not logged (ana only locates the Ohmic cell). Candidates: Ohmic heating or
  field diffusion into the floored top (rho = 1e-13) raising c_f there. Worth one diagnostic bin at the dt-min cell.
- b3_e11 is the costliest 3 G no-Ohmic arm late (11.6 min/rot vs 8.4 for e12): same cause as above, and it is a
  single realisation (twin noise on late dt: 3.02 vs 3.12 s = 3 %).
- eta_num(est) is negative in every arm and scales with the cap (-2e10 at e0 to -1e12 at e13): the missing
  boundary Poynting flux / grid-scale J dominates; not usable as a numerical-diffusivity estimate (also local
  change: dME/dt window 0.25 rot, `ana_eta.py.orig` is the shipped one).
- Memory said b10_e13_sts stopped at "99.6 %"; the log gives t = 3.325791e7 = 96.1 % (rot 301.92). The waiter's
  "~9.4k cycles" matches 96.1 %.
- rotation.txt has no b3_e0_f13 line (first link list is in rotation.old; final job 3272979 per done.txt).
  queue.txt still lists 4 ranks per arm (pre 2-GPU switch); the logs show 2 (twin 1).
- The viper arms use binary ae767d20, which predates the seam-flux request fix 2dfe4e95. The leak aborted Cray
  MPICH at cycle 8100 on 2 ranks here; if viper's MPI does not abort it may still leak requests over the ~1e5-1e6
  cycles of the e14 / f16 arms. Check before releasing 12028523-25.
- div B: b3_e13_sts 2.7e-11 is ~3x the explicit twin (1.0e-11); both round-off level, but STS adds stages.

## Files
- Per arm: `runs/<arm>/ana_eta.txt` (full band tables), `runs/<arm>/run.<jobid>.log`, `logs/w<arm>.<jobid>.out`.
- `collect_arm.sh ARM [OHMIC_DT]`: stability counts, per-log dt/wall/ms-per-cycle, last-0.5-rot rates, RKG stage
  histogram, TOTAL over restart links, and the ana_eta summary lines.

## Appendix: ana_eta.txt summary lines (verbatim)

```
# arm runs/b3_e0_f13: bbot 8.462844e-01 (code), max_eta 1.0e+00, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 6.5 min/rot, 14.46 ms/cycle, dt mean 4.116 s min 2.913 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 9.639 s, Ohmic estimate 17462333414163.836 s at p 7.49 bar lat -35 lon +136 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 3.765e+07 erg/s = 1.24e-23 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 1.761 s, Ohmic estimate 17462333414163.836 s at p 7.53 bar lat -35 lon +136 -> limiter CFL/other
   jet max 13073 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8576/5071/2351 m/s; T(0.1 bar) day 2048 night 1482 contrast 566 K
   E_mag 2.356e+33 E_kin 1.217e+35 erg, |B|max 363; Q_ohm 7.422e+17 erg/s = 2.44e-13 L_abs; W_lorentz -7.971e+27; dEM/dt 1.252e+28; eta_num(est) -2.76e+10 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 6.32e-12
```

```
# arm runs/b3_e11_f13: bbot 8.462844e-01 (code), max_eta 1.0e+11, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 7.7 min/rot, 14.32 ms/cycle, dt mean 3.414 s min 1.960 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 9.639 s, Ohmic estimate 174.623 s at p 7.49 bar lat -35 lon +136 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 1.441e+18 erg/s = 4.73e-13 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 2.020 s, Ohmic estimate 174.623 s at p 7.53 bar lat -35 lon +136 -> limiter CFL/other
   jet max 12767 m/s at 0.00075 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8757/5046/2353 m/s; T(0.1 bar) day 2048 night 1479 contrast 569 K
   E_mag 1.718e+33 E_kin 1.221e+35 erg, |B|max 194; Q_ohm 9.253e+27 erg/s = 3.04e-03 L_abs; W_lorentz -1.210e+28; dEM/dt 8.695e+27; eta_num(est) -2.86e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 7.02e-12
```

```
# arm runs/b3_e12_f13: bbot 8.462844e-01 (code), max_eta 1.0e+12, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 7.4 min/rot, 14.30 ms/cycle, dt mean 3.562 s min 2.402 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 9.639 s, Ohmic estimate 17.462 s at p 3.91 bar lat +48 lon +165 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 1.074e+19 erg/s = 3.53e-12 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 3.263 s, Ohmic estimate 17.462 s at p 3.99 bar lat +48 lon +162 -> limiter CFL/other
   jet max 12881 m/s at 0.0013 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9216/5066/2354 m/s; T(0.1 bar) day 2044 night 1480 contrast 564 K
   E_mag 1.028e+33 E_kin 1.229e+35 erg, |B|max 118; Q_ohm 1.003e+28 erg/s = 3.29e-03 L_abs; W_lorentz -1.181e+28; dEM/dt 4.136e+27; eta_num(est) -7.32e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 7.56e-12
```

```
# arm runs/b3_e12_f13_r2: bbot 8.462844e-01 (code), max_eta 1.0e+12, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 9.8 min/rot, 18.97 ms/cycle, dt mean 3.569 s min 2.256 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 9.639 s, Ohmic estimate 17.462 s at p 3.91 bar lat +48 lon +165 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 1.074e+19 erg/s = 3.53e-12 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 3.061 s, Ohmic estimate 17.462 s at p 3.99 bar lat +48 lon +162 -> limiter CFL/other
   jet max 12763 m/s at 0.0013 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9135/5066/2354 m/s; T(0.1 bar) day 2044 night 1480 contrast 564 K
   E_mag 1.023e+33 E_kin 1.229e+35 erg, |B|max 118; Q_ohm 1.004e+28 erg/s = 3.30e-03 L_abs; W_lorentz -1.200e+28; dEM/dt 4.099e+27; eta_num(est) -7.37e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 4.85e-12
```

```
# arm runs/b3_e13_f13: bbot 8.462844e-01 (code), max_eta 1.0e+13, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 14.7 min/rot, 14.13 ms/cycle, dt mean 1.760 s min 1.755 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 1.778 s, Ohmic estimate 1.778 s at p 1.8 bar lat +52 lon +171 -> limiter Ohmic
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 8.636e+19 erg/s = 2.83e-11 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 0.875 s, Ohmic estimate 1.755 s at p 1.67 bar lat +59 lon +163 -> limiter CFL/other
   jet max 13606 m/s at 0.00075 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9576/5089/2362 m/s; T(0.1 bar) day 2041 night 1477 contrast 564 K
   E_mag 3.808e+32 E_kin 1.238e+35 erg, |B|max 101; Q_ohm 7.733e+27 erg/s = 2.54e-03 L_abs; W_lorentz -9.754e+27; dEM/dt 8.265e+26; eta_num(est) -1.03e+12 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 1.04e-11
```

```
# arm runs/b3_e13_sts_f13: bbot 8.462844e-01 (code), max_eta 1.0e+13, use_rkg_sts true, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 17.1 min/rot, 25.55 ms/cycle, dt mean 2.742 s min 2.066 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 9.639 s, Ohmic estimate 1.778 s at p 1.8 bar lat +52 lon +171 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 8.636e+19 erg/s = 2.83e-11 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 2.560 s, Ohmic estimate 1.755 s at p 1.67 bar lat +59 lon +163 -> limiter CFL/other
   jet max 13692 m/s at 0.00075 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9585/5088/2362 m/s; T(0.1 bar) day 2041 night 1477 contrast 564 K
   E_mag 3.801e+32 E_kin 1.238e+35 erg, |B|max 101; Q_ohm 7.547e+27 erg/s = 2.48e-03 L_abs; W_lorentz -9.701e+27; dEM/dt 8.802e+26; eta_num(est) -1.02e+12 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 2.74e-11
```

```
# arm runs/b10_e0_f13: bbot 2.820948e+00 (code), max_eta 1.0e+00, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 14.9 min/rot, 14.09 ms/cycle, dt mean 1.731 s min 1.055 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 3.090 s, Ohmic estimate 17462333414163.836 s at p 7.49 bar lat -35 lon +136 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.987e+30 E_kin 1.278e+35 erg, |B|max 2.8; Q_ohm 4.184e+08 erg/s = 1.37e-22 L_abs; W_lorentz -1.858e+23; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 0.575 s, Ohmic estimate 17462333414163.836 s at p 7.52 bar lat -35 lon +136 -> limiter CFL/other
   jet max 11395 m/s at 0.00075 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8079/4923/2314 m/s; T(0.1 bar) day 2064 night 1477 contrast 588 K
   E_mag 7.827e+33 E_kin 1.173e+35 erg, |B|max 520; Q_ohm 3.060e+18 erg/s = 1.00e-12 L_abs; W_lorentz -2.334e+28; dEM/dt 3.632e+28; eta_num(est) -1.95e+10 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 9.28e-12
```

```
# arm runs/b10_e11_f13: bbot 2.820948e+00 (code), max_eta 1.0e+11, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 16.8 min/rot, 14.20 ms/cycle, dt mean 1.551 s min 0.936 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 3.090 s, Ohmic estimate 174.623 s at p 7.49 bar lat -35 lon +136 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.987e+30 E_kin 1.278e+35 erg, |B|max 2.8; Q_ohm 1.601e+19 erg/s = 5.26e-12 L_abs; W_lorentz -1.858e+23; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 0.570 s, Ohmic estimate 174.623 s at p 7.52 bar lat -35 lon +136 -> limiter CFL/other
   jet max 11398 m/s at 0.001 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8167/4946/2325 m/s; T(0.1 bar) day 2067 night 1483 contrast 584 K
   E_mag 5.170e+33 E_kin 1.182e+35 erg, |B|max 304; Q_ohm 2.418e+28 erg/s = 7.93e-03 L_abs; W_lorentz -3.446e+28; dEM/dt 2.221e+28; eta_num(est) -2.88e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 9.54e-12
```

```
# arm runs/b10_e12_f13: bbot 2.820948e+00 (code), max_eta 1.0e+12, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 18.7 min/rot, 14.13 ms/cycle, dt mean 1.384 s min 0.911 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 3.090 s, Ohmic estimate 17.462 s at p 3.91 bar lat +48 lon +165 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.987e+30 E_kin 1.278e+35 erg, |B|max 2.8; Q_ohm 1.194e+20 erg/s = 3.92e-11 L_abs; W_lorentz -1.858e+23; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 0.920 s, Ohmic estimate 17.462 s at p 3.99 bar lat +48 lon +162 -> limiter CFL/other
   jet max 11197 m/s at 0.0042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8065/5016/2334 m/s; T(0.1 bar) day 2057 night 1486 contrast 571 K
   E_mag 3.216e+33 E_kin 1.198e+35 erg, |B|max 196; Q_ohm 2.490e+28 erg/s = 8.17e-03 L_abs; W_lorentz -3.222e+28; dEM/dt 1.119e+28; eta_num(est) -6.48e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 6.57e-12
```

```
# arm runs/b10_e13_f13: bbot 2.820948e+00 (code), max_eta 1.0e+13, use_rkg_sts false, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 23.8 min/rot, 14.18 ms/cycle, dt mean 1.095 s min 0.848 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 1.778 s, Ohmic estimate 1.778 s at p 1.8 bar lat +52 lon +171 -> limiter Ohmic
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.987e+30 E_kin 1.278e+35 erg, |B|max 2.8; Q_ohm 9.596e+20 erg/s = 3.15e-10 L_abs; W_lorentz -1.858e+23; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 0.752 s, Ohmic estimate 1.755 s at p 1.67 bar lat +59 lon +163 -> limiter CFL/other
   jet max 11676 m/s at 0.0032 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8785/5092/2350 m/s; T(0.1 bar) day 2047 night 1486 contrast 561 K
   E_mag 1.486e+33 E_kin 1.215e+35 erg, |B|max 192; Q_ohm 2.040e+28 erg/s = 6.69e-03 L_abs; W_lorentz -2.447e+28; dEM/dt 1.476e+27; eta_num(est) -7.45e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 1.52e-11
```

```
# arm runs/b10_e13_sts_f13: bbot 2.820948e+00 (code), max_eta 1.0e+13, use_rkg_sts true, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 42.0 min/rot, 25.44 ms/cycle, dt mean 1.113 s min 0.827 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 3.090 s, Ohmic estimate 1.778 s at p 1.8 bar lat +52 lon +171 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.987e+30 E_kin 1.278e+35 erg, |B|max 2.8; Q_ohm 9.596e+20 erg/s = 3.15e-10 L_abs; W_lorentz -1.858e+23; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00159.bin rot 302.000: dt run 0.827 s, Ohmic estimate 1.755 s at p 1.67 bar lat +59 lon +163 -> limiter CFL/other
   jet max 11607 m/s at 0.0032 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 8668/5094/2350 m/s; T(0.1 bar) day 2047 night 1486 contrast 562 K
   E_mag 1.490e+33 E_kin 1.214e+35 erg, |B|max 190; Q_ohm 2.046e+28 erg/s = 6.72e-03 L_abs; W_lorentz -2.486e+28; dEM/dt 1.549e+27; eta_num(est) -7.52e+11 cm2/s
# div B (dhj.mhd_divb.00003.bin): max |divB| dr/|B| = 1.07e-10
```

```
# arm runs/b3_e14_sts_f13: bbot 8.462844e-01 (code), max_eta 1.0e+14, use_rkg_sts true, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 32.7 min/rot, 42.16 ms/cycle, dt mean 2.370 s min 1.873 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 9.639 s, Ohmic estimate 0.185 s at p 0.89 bar lat +52 lon +171 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.188e+29 E_kin 1.278e+35 erg, |B|max 0.841; Q_ohm 7.232e+20 erg/s = 2.37e-10 L_abs; W_lorentz -1.672e+22; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00158.bin rot 302.000: dt run 1.946 s, Ohmic estimate 0.181 s at p 0.893 bar lat +58 lon +178 -> limiter CFL/other
   jet max 13819 m/s at 0.00075 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9584/5086/2363 m/s; T(0.1 bar) day 2041 night 1477 contrast 564 K
   E_mag 2.482e+32 E_kin 1.239e+35 erg, |B|max 95.8; Q_ohm 5.732e+27 erg/s = 1.88e-03 L_abs; W_lorentz -8.353e+27; dEM/dt 2.641e+26; eta_num(est) -9.39e+11 cm2/s
# div B (dhj.mhd_divb.00002.bin): max |divB| dr/|B| = 1.63e-10
```

```
# arm runs/b10_e14_sts_f13: bbot 2.820948e+00 (code), max_eta 1.0e+14, use_rkg_sts true, rsolver lhlld
# cost: 2.00 rot run (to rot 302.00), wall 60.9 min/rot, 32.21 ms/cycle, dt mean 0.971 s min 0.626 s
## dhj.mhd_w_bcc.00151.bin rot 300.000: dt run 3.090 s, Ohmic estimate 0.185 s at p 0.89 bar lat +52 lon +171 -> limiter CFL/other
   jet max 14207 m/s at 0.00042 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9775/5075/2331 m/s; T(0.1 bar) day 2039 night 1481 contrast 558 K
   E_mag 7.987e+30 E_kin 1.278e+35 erg, |B|max 2.8; Q_ohm 8.035e+21 erg/s = 2.64e-09 L_abs; W_lorentz -1.858e+23; dEM/dt nan; eta_num(est) nan cm2/s
## dhj.mhd_w_bcc.00160.bin rot 302.000: dt run 0.675 s, Ohmic estimate 0.181 s at p 0.893 bar lat +58 lon +178 -> limiter CFL/other
   jet max 12785 m/s at 0.001 bar; u_eq(|lat|<20) 1e-3/0.1/1 bar = 9303/5113/2354 m/s; T(0.1 bar) day 2045 night 1485 contrast 560 K
   E_mag 1.122e+33 E_kin 1.222e+35 erg, |B|max 155; Q_ohm 1.663e+28 erg/s = 5.46e-03 L_abs; W_lorentz -2.156e+28; dEM/dt 1.158e+27; eta_num(est) -7.16e+11 cm2/s
# div B (dhj.mhd_divb.00004.bin): max |divB| dr/|B| = 2.01e-10
```

```
# arm runs/b10_e14_f13: TBD(b10_e14_f13) (viper arm, not run on DeltaAI)
```

```
# arm runs/b10_e13_f16: TBD(b10_e13_f16) (viper arm, not run on DeltaAI)
```

```
# arm runs/b3_e13_f16: TBD(b3_e13_f16) (viper arm, not run on DeltaAI)
```
