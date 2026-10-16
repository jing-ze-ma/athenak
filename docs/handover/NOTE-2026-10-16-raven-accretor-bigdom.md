# NOTE 2026-10-16 raven -> caltech: accretor bigdom v3 (pinned-rate nozzle at 14.5), 2-D test - INTERIM

Reply to TASK-2026-10-15-raven-accretor-bigdom.md **v3 (1db9e3e7)**. The static phase and the stream phase up to
t 0.02506 have run. That is just before the first impact (the stream front is at r 10.3 at t 0.0251). The rest of the
run, to tlim 0.036 (impact, dt retries, arrival phi), is queued on the Raven gpu partition (job 31212975). It
continues from the newest stream rst. Raven work dir: /raven/ptmp/jinma/accbig_1016 (v3/ = this test).

## Build and runs
- Code: accretor-bigdom-1015 @ efe33e6a (code 7fb237a8). Worktree /raven/u/jinma/ATHENAK/wt-accbig.
- Build: gcc/13 + cuda/12.6 + openmpi_gpu/5.0, Kokkos 4.6.02 CUDA AMPERE80, Release, MPI, `-D PROBLEM=ry_per_accretor
  -D CMAKE_CXX_FLAGS=-ffp-contract=off`, full fresh build (job 31212898).
  Binary `bin/athena_a100_ry_per_accretor_efe33e6a_nofma`, **md5 b436ea574252f57166e9f24be575df88**.
  nvcc_wrapper passes -ffp-contract=off to the host compiler only; device code keeps nvcc's default --fmad=true.
  This is the same as the Caltech nofma recipe.
- Deck: v3 bigdom.athinput (md5 3a818e64) with @HDIR@ replaced and `dt_report = true` added. Tables md5 OK.
- Layout: 4 A100 (gpudev, 1 node), 4 ranks x 2 MeshBlocks (840x4x256). The TASK said 1 GPU. On the Raven gpu
  partition a 1-GPU job would wait ~1.5 d, and 1 A100 costs 1.8-2.0 s/cycle (v2 static, 418 cycles), so 4 A100 it is.
- Jobs: static, gpudev link 1 = 31213087 (t 0 -> 0.00946, 1356 cycles). Stream, gpudev link 2 = 31213158
  (0.00946 -> 0.01905) and link 3 = 31213159 (-> 0.02506, cycle 4795). Remainder: gpu partition 31212975
  (4 A100, 3 h, `jobs/chain.sbatch`, starts from the newest v3/stream rst; Slurm estimate ~8 h out at writing).
- Side checks: 31212974 (v3 vs v2 static bitwise, 50 cycles, 1 A100); 31213160 (orbit1 fdb18517 static to t 0.0012
  for the envelope comparison); v2 regression gate ddb14cf8 vs fdb18517 BITWISE (31212550/31212563, s5 deck, 20 cycles).
  The v2 static run was not reused: v3 static differs from v2 static from cycle 1 (Eirr, Mmsk, Ltop and the 1-mom
  hst columns at 1e-8). dt is identical, but the bins and hst are not bitwise equal. So the static phase started fresh on v3.

## Results vs pass checks (t <= 0.02506; the impact itself is pending)
| check | result | verdict |
|---|---|---|
| No FATAL / NON-CONVERGED storms; dt_retry list | No FATAL. The static phase has 0 NC and 0 retries. **8 dt retries, all just after the switch-on**: cycles 1366-1431, t 0.009516-0.009585, 7 x eint<=0 and 1 Picard NC pass (resid 4.9e15, cycle 1388). All recovered by cycle 1448. The bad cells are r 14.15-14.19, phi -0.1..1.3 deg: the fresh stream front 0.35 Rsun below the nozzle hitting the floor gas. No retries after that (link 3: 0). Picard mean 2.1-2.3, max 5 (53 once at switch-on). | PASS (switch-on transient, recovered) |
| Pinned rate, Mnzf (col 25) = 19.795 within 0.5 % | 19.7989 at the first output after switch-on (t 0.0105), falling slowly to 19.7969 at t 0.019 (**+0.01 to +0.02 %**). | PASS |
| Min - Mnzf constant after the fill | 0.016209 -> 0.016223 over t 0.0105-0.0190: drift 1.7e-3 per time unit (1e-4 of the rate). | PASS |
| Nozzle stable, no growing junk at r 14.4-14.7, \|phi-2\| < 25 deg | Nozzle rho max 0.2726 (= the pinned peak), the ridge at 14.4 is phi 2.2 deg, and the stream leaves cleanly. No gas above floor in 14.7-16.3 or in the lobe. **One artefact**: the first back-sink row behind the nozzle (r 14.607, phi -12..+16 deg) holds floor density with **v_r = +2e5 to +9.7e5 km/s** in every hydro_w bin (falling with time: 8.5e5 at t 0.0106, 4-5e5 at 0.025). The sink resets these cells to floor at rest every stage (ry_per_accretor.cpp:1150-1160), and the hydro dt does not see them (dt from bins would be 1e-8). So this is momentum the floor gas gains after the reset within the step (radiation feedback or the radeq keep-momentum floor). It is harmless to dt and mass, but it shows in outputs and Emsk. Recommend: zero v (or skip feedback) in sink cells after the radiation step, or exclude the row from the output. | PASS with artefact |
| Mmsk small, smooth | Static 2.5e-10 by t 0.00946. Stream link 2: 1.7e-5 by t 0.019, smooth, rate falling (1.6e-3 -> 1.05e-3 per time unit). | PASS |
| Same stream path as orbit1 (arrival phi 15-16.5 deg) | Ridge at t 0.0251: r 14.4 phi 2.2, 13.8 phi 3.3-3.6 (orbit1 window 3.41), 13.0 4.7, 12.0 6.6, 11.0 8.2, 10.3 9.2 deg (front). The first-impact location and time are **PENDING** (gpu job). The front speed is consistent with the TASK estimate (13.8 at t ~0.016 = switch-on + 0.0065; impact ~0.031). | PENDING (path so far consistent) |
| Mask stable, lobe at floor, donor light sensible | Lobe (21438 cells) is at rho 1e-8 with \|v\| < 0.03 km/s in every bin, and nothing hot leaks (r > 16.3 outside the lobe: rho <= 1.17e-8, \|v\| < 2 km/s). Eirr rate: static 7.37e6 code/time, against sigma T_d^4 x exposed lobe area in the band (7.23 interface + 5.13 r_out face = 12.4 Rsun^2) x 5.09e5 = 6.29e6, so 1.17x: **sensible**. Stream phase: 1.55e7, because the back-sink sector (r > 14.59, \|dphi\| < 4 sigma) also holds E = a T_d^4 (inmk includes it, ry_per_accretor.cpp:1135), so the region between the nozzle and L1 shines at 26 kK too: a design point to decide. Ltop 2.0e5 static (orbit1 4.5e4), 3.5-4.3e5 in the stream phase. Note: at output time the lobe m1 E is 1.1-3.4 (16-50 % of a T_d^4 = 6.787) with \|F\| ~ cE: outputs come after the implicit radiation step that follows the reset. | PASS (2 design notes) |
| dt not set by the nozzle or the mask; reference 7.03-7.05e-6 | **dt is not set by the nozzle or the mask, but it does drop: 7.03e-6 until t 0.008, 6.43e-6 at the end of the static phase, then 5.9e-6 (t 0.0105), 7.03e-6 (t 0.0125-0.015), 5.5e-6 (0.017), 4.3e-6 (0.020), 2.8e-6 (0.0225), 2.4-3.3e-6 (t 0.0245-0.025).** It is set by thin atmosphere gas above the star: rho 3e-7..3e-6, r 9.1-9.6, at phi ~87-97 deg and ~263-268 deg (both sides, 90 deg from the donor axis), with v_phi +700 / -940 km/s (both flowing away from the donor side toward phi 180) and v_r up to +180 km/s. Its peak speed grows steadily: ~200 km/s at t 0.0095 (static, phi -50..-60 deg), 370 at 0.015, 630 at 0.020, 940 at 0.025 (local escape speed ~800 km/s). orbit1 (no donor light) had 7.03e-6 to impact, so **the most likely driver is donor_irr light pushing the thin atmosphere around the star**. It is unconfirmed: the clean test is a static arm with donor_irr = false (a 4xA100 gpudev link, ~12 min); not run (link budget). Caveat for dt_report: it prints rank 0's local hydro minimum only (ry_per_accretor.cpp:1367, my_rank == 0). On >1 rank the global cell can be elsewhere, so I located the cell from the bins (approximate gas CFL, tool/dtcell.py). | **FAIL-ish / FLAG**: not nozzle/mask, but dt 2-3x below the reference before impact |
| Envelope r <= 10.4 unchanged vs orbit1 (static) | Same face grid to 0 (r <= 10.4018). At t 0.00118 (cycle 168, both runs): r <= 8.6 agrees to 3e-6 (rho), gas with rho > 0.01 agrees to median 5e-6 relative (max 0.15 at the photosphere edge). The thin atmosphere r 9.0-9.6 (rho < 1e-2) differs at O(1): both runs have fast thin flows there (orbit1: 84, -140 km/s at rho 1e-7), so the differences come from the new outer radiation (donor light, dark top). | PASS for the envelope proper; the atmosphere differs (see dt) |
| Cost per cycle | 4 A100: static 0.458 s/cycle, stream 0.41-0.45 s/cycle. 1 A100 (v2 static, same grid): 1.8-2.0 s/cycle, i.e. ~2.7x the orbit1 H200 0.51-0.56 x 1.31 cells. | - |
| Time of first impact | PENDING (gpu job 31212975). | PENDING |

## Recommendations
1. Test donor_irr = false (static, ~1350 cycles) to confirm or refute that the donor light drives the 900 km/s
   thin-atmosphere flow that halves dt before impact. If confirmed: limit the radiation force on floor/atmosphere
   gas, or move the light source.
2. Back-sink row artefact (r 14.607): zero the momentum of sink cells after the radiation/floor step.
3. Decide whether the back-sink sector should radiate at T_d (it doubles Eirr and lights the stream from behind).
4. dt_report: reduce over ranks (MPI MINLOC) so that multi-rank runs report the global cell.

Scripts: /raven/ptmp/jinma/accbig_1016/jobs/{chain,bw,env}.sbatch, tool/{ana_v3,dtcell,envloc}.py.

## ADDENDUM 1 (2026-10-16 morning): donor_irr = false diagnostic - the donor light drives the dt drop (CONFIRMED)
Job 31214095 (gpudev, 4 A100, 638 s): the v3 deck with `problem/donor_irr=false`, fresh static phase to t 0.00946 (1345
cycles), same binary efe33e6a. Dir /raven/ptmp/jinma/accbig_1016/v3/noirr. Compared with the v3 static run (donor_irr on).
- **dt**: no irradiation: 7.0473e-6 -> 7.0291e-6 flat to the end (0 retries, 0 NC, Picard mean 2.02). With irradiation:
  6.90e-6 at cycle 1200 and 6.43e-6 at cycle 1300. **The static dt drop is caused by donor_irr.**
- **Thin atmosphere r 9.1-9.6** (equatorial pair, rho < 1e-4 for the mean):
  | t | max\|v\| irr | max\|v\| no irr | mean\|v_phi\| irr | no irr |
  |---|---|---|---|---|
  | 0.00237 | 67 | 12 | 6.5 | 1.1 |
  | 0.00473 | 116 | 22 | 23.6 | 1.9 |
  | 0.00710 | 172 | 40 | 58.7 | 3.1 |
  | 0.00946 | 228 | 47 | 62.2 | 5.6 |
  This is 5x the speed. orbit1-like slow atmosphere motions remain without irradiation (max 47 km/s); 200-940 km/s needs
  the donor light.
- **How much donor light reaches r 9-10 (m1 bins at t 0.00946):**
  - E_rad rises by a nearly uniform **factor 4.7-5.6 at all radii 9.2-16 and at all phi, including the far side
    (phi 180: 4.85x)**. Median E at r 9.6-10.4 is 13.7 with irradiation vs 2.8 without, i.e. **2x the donor's own a T_d^4
    (6.787)** and above donor + star (9.6). At r 9.0-9.6 it is 17.8 vs 3.5.
  - \|F\| ~ cE (4.6e6 vs cE 4.1e6, free streaming), mostly AZIMUTHAL: dF_phi = +3.6..+4.5e6 at phi 30-90 and 210-240,
    -3.4..-4.5e6 at phi 120-150 and 270-330 (r 9.6). That is 5-20x the star's own radial flux there (median 2.2e5).
  - Expected donor flux on the star's facing side: real 3-D dilution L_d/(4 pi D^2), D = a - R ~ 23.8, gives 1.5e5 code. An
    in-band view-factor estimate (lobe subtends ~40 deg) gives ~1-2e5. The model's excess flux at r 9.6 is **~30x larger
    and present on the night side too**, so the donor irradiation as implemented over-irradiates the star by over an
    order of magnitude.
  - Budget: injected Eirr 6.87e4 by t 0.00946 (7.3e6 per time unit). Excess radiation energy resident in the domain
    is only ~1.0e3 (r 9-24), so it is not a trapped reservoir. It is a strong, continuously replenished azimuthal flux.
  - **Back-sink sector**: NOT needed for the effect. The whole effect is present in the static phase, where the back
    sink is off (it needs the stream on). In the stream phase the back sector also shines at T_d and doubles Eirr
    (1.55e7), and the atmosphere speeds kept growing (630 -> 940 km/s), so it adds to the problem.
- **Probable causes (no code changed):**
  (a) band geometry: the +-0.35 deg band with reflecting theta walls has no vertical dilution (in-plane 1/d instead of
  1/d^2), so a lobe held at a T_d^4 over its full in-band surface is far brighter at the star than a 3-D donor would be;
  (b) closure = vet_col builds the Eddington tensor from RADIAL columns. A lateral source (the donor light comes in from
  the side) cannot be represented: the uniform ~5x E at all phi and the strong azimuthal F (> star's flux) point to the
  closure, not to geometry alone.
- **Proposed fixes (for Caltech to choose; not implemented):**
  1. Do not irradiate from the back-sink sector (mask only).
  2. Replace "lobe holds E = a T_d^4" with a diluted, directed source: an inward flux of the 3-D value L_d/(4 pi d^2)
     (d from the donor centre) entering through the r_out boundary (or the lobe surface) only within the donor's solid
     angle. Or keep the lobe but scale T_d so that the band flux at the star's sub-donor point equals the 3-D value
     (here ~1/30 of the present flux, i.e. T_d,eff ~ 26 kK x 30^-1/4 ~ 11 kK).
  3. Check the vet_col closure with a lateral beam (a short test with closure = m1 or the VET recomputed from the full
     angular field) before trusting any irradiated run.
  4. As a guard: limit the radiation force on floor/atmosphere gas (rho < ~1e-4 code) or extend the floor sponge to
     the atmosphere above the photosphere.
  The v3 gpu job (31212975) is left running unchanged as instructed: it tests v3 as specified, impact included.

Ack d7b8e051 (link STOP rule): Raven bsgnc_raven.sh / e1chain_raven.sh and Viper link_viper.sh / e1chain_viper.sh already stop only on the summary count NON-CONVERGED=[1-9] (+FATAL/DIVERGED/NaN/rc); DT_RETRY reasons never contain DIVERGED; no change needed; rule recorded.
(The accbig chain script jobs/chain.sbatch of 31212975 does not test NON-CONVERGED at all: it continues a phase only after "Terminating on time limit".)

## ADDENDUM 2 (2026-10-16 evening): v3 through impact (gpu job 31212975, 4 A100, 53 min) - over-irradiated v3 deck, for comparison
Restart from the link-3 rst (t 0.02506), the v3 deck unchanged (donor_irr on at full strength), tlim 0.036 reached
(cycle 13258, 8463 cycles in 3189 s = 0.37 s/cycle).
- **No FATAL; 0 NON-CONVERGED; 0 dt retries in this link** (the whole v3 run has 8 retries, all at the switch-on).
  orbit1 needed ROL2U + 21 retries at its impact; here dt was already 2-6x smaller when the stream arrived.
- **dt through impact**: 2.8-3.3e-6 (t 0.025), 1.53e-6 from t 0.026 (stream infall onto the atmosphere, v_r ~ -500
  km/s at rho 1e-5..1e-2), 1.15e-6 from t 0.0315, 0.84-1.03e-6 at t 0.034-0.036. The rank-0 dt_report cells are
  stream/atmosphere gas at r 9.4-9.6 (v_r -330..-510, v_phi 200-680 km/s).
- **First impact**: dense infalling stream gas (v_r < -50, rho > 1e-3) first reaches r 9.4-9.65 at **t ~0.027** (587
  cells at t 0.0272, phi 0.6-18.9 deg; 8010 at t 0.0284). Shock heating of dense gas starts at **t ~0.031-0.032**:
  hottest dense gas at r 9.30 with e/rho 3294 (t 0.0319) -> 4044 (t 0.036; the far side stays at 2470), at **phi 14.0 ->
  15.2 deg**, v_r -346..-378 km/s.
- **Arrival phi**: density ridge at r 9.4-9.6 at phi 13.4-14.3 deg (t 0.033-0.036); hottest impact gas at phi 14.0-15.2
  deg. **Within/at the lower edge of the orbit1 window (15-16.5 deg; orbit1 hottest gas r 9.43, phi 15 deg).** The
  impact time (~0.031 for the heating) matches the TASK estimate.
- This deck over-irradiates (ADDENDUM 1), so the dt numbers are not representative. v4 (scaled, gravity-darkened
  donor light) is being tested now (NOTE-2026-10-16/17-raven-accretor-v4.md).
