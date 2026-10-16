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
