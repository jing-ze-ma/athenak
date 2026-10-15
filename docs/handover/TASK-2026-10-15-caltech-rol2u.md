# TASK raven -> Caltech: ROL2U from the newest production rst of the He giant (hgF4w) and AG Car E1 (agcE1s) — gate, short BASE window, production-capable ROL2U chains

From Raven worker CR2 (/raven/ptmp/jinma/caltech_rol2u_1015) for the coordinator to push. Bundle: `docs/handover/rol2u-1015/caltech_rol2u/` (copy next to
this TASK into docs/handover/caltech_rol2u-1015/; on Caltech put it at `/resnick/groups/carnegie_poc/jingze/rol2u_1015/caltech_rol2u`).

USER (10-15): "try ROL2U on Caltech on the latest restart files of the He giant and AG Car". USER update (via coordinator, 10-15): the
ROL2U arms are PRODUCTION-CAPABLE chains, not short tests (full production-style chain from the newest prod rst, tlim 1e9, 7-d links,
production hooks/guards/cadence, own run dir), so a ROL2U chain can simply BECOME the production run if it runs fine (user decision
then). A BASE arm runs only as a short window (same start) for speed/accuracy and as the FIN4-rst readability / bitwise gate.
ROL2U must pass the same gate (smoke + restart bitwise) before its long links start.

**Do NOT touch production**: hgF4w (hegiant_1014/wide4/run; queued link 4513229/30, guard 4514966/67) and agcE1s (e1s_1014/prod,
4513330-32) keep running/queued unchanged. Production dirs are READ only (the rst is COPIED once by seed_cr2.sh; never hard-link,
move or write anything there). No production switch: that is a later user decision.

## What ROL2U is (fin51-r @ 1bd0c2ce, origin/fin51-r)
`fin51-r` = FIN5.1 (fin51-1014 @ 21c4839d, your validated build) + the R keys (29e0d699) + OL2U keys (1bd0c2ce); every new key is
default-off: 1bd0c2ce with default keys is BITWISE to 21c4839d (Raven A100 E1 30 cyc rst + 31 hst; CPU IDORT 20/20).
ROL2U = these keys under `<rad_m1>` (exact names, verified in src/rad_m1/rad_m1_implicit.cpp:1000-1001, 1810-1824):
```
rad_m1/implicit_tol_rel=0.1  rad_m1/implicit_tol_rel_max=1e-4  rad_m1/implicit_opac_lag=pred  rad_m1/implicit_opac_lag_passes=2
rad_m1/implicit_gas_newton_u4=true
```
Validation on Raven: full IDORT (292 dirs x 3 arms) PASS all families (ridort_1014 NOTE); long drift on AG Car P268/S3 (8e6 s each,
A100/MI300A): ROL2U - BASE inside the cfl/2 and chaos yardsticks, Picard 29 -> 5.3, s/cycle -55..-69 %, 0 DT_RETRY; R-only effect:
more solved-E <= e_floor clips (S3 2009 vs 140 cell-solves; energy 4e-6 of L0 t). The He giant has NO ROL2U run yet anywhere (Raven/Viper
H_* arms still queued): the gate below is its first test.

## 0. Restart compatibility FIN4 / fin4-dtretry rst -> 1bd0c2ce (checked from code + Raven pre-test)
- rst FORMAT unchanged from FIN4 (d505d21f) and 5b7e1bc5 to 1bd0c2ce: the FIN5 "chunked rotavg rst pack" (RotAvgRstPackHost) writes the
  same bytes; the reader changes are int -> size_t casts only (pgen.cpp tail offsets, restart.cpp data_size); the hesdirk2-slope /
  predictor / rotavg slabs are read as before. dt_retry adds no binary trailer, only `<time>` keys (dt_retry_cap, dt_retry_ngood_now,
  written only while a cap is or was set) -> 5b7e1bc5 E1 rst carry them and FIN5.1 restores them (Raven/Viper G = G2 BITWISE).
- dt_retry default on a RESTART is OFF unless the key is in the rst/deck/command line (dt_retry.cpp:185-189). E1 rst: key present
  (true). He giant FIN4 rst: key ABSENT -> give `time/dt_retry=true` on the command line: on a restart an absent NAME in an existing
  block is ADDED (parameter_input.cpp ModifyFromCmdline add_absent; the "absent key FATAL" rule applies to FRESH starts only).
- Restart header limit (your 6fee27be finding): ParameterInput reads the rst header in 4096-byte chunks and FATALs if `<par_end>` is not
  within 45 056 bytes. Your 200-Rsun He giant FIN4 rst header is 41.6 KB; FIN5.1 dumps ~0.4 KB more, ROL2U + dt_retry keys ~0.5 KB more
  -> ~42.5 KB expected: fits, but the gate CHECKS the <par_end> offset of the ROL2U rst (must stay < 44 000 to leave margin). E1S 38.9 KB.
- Expected: BASE (1bd0c2ce, production keys, from the FIN4 / 5b7e1bc5 rst) is BITWISE to the production binary continuing from the same
  rst (FIN5 = FIN4 and FIN5.1 = FIN5 / = 5b7e1bc5 BITWISE on A100 + MI300A + your H200 smoke; the FIN4 He giant fused-twin IC overflow
  affects FRESH starts only). This is GATE-BASE below.
- Raven pre-test (A100, 1 node x 4, same code; FIN4 5c504e49 / 5b7e1bc5-written rst restarted by both binaries; jobs/xrst.sbatch,
  jobs/xrol.sbatch in /raven/ptmp/jinma/caltech_rol2u_1015): RESULTS_PLACEHOLDER

## 1. Build fin51-r @ 1bd0c2ce (H200, nofma) — record md5 + `git rev-parse HEAD`
In your FIN5.1 tree (`/resnick/groups/carnegie_poc/jingze/builds_fin51`, built 21c4839d): `git fetch origin fin51-r && git checkout 1bd0c2ce`
(4 files differ from 21c4839d: rad_m1.hpp, rad_m1_dtretry_visit.hpp, rad_m1_implicit.cpp, rad_m1_opacity.hpp; no new source file ->
incremental build). EXACTLY the fin51 recipe (bsg_1009/build_inc_nofma.sh with BUILD_ROOT=builds_fin51): gcc 13.2.0 + cuda 12.9.0 +
hpcx 2.17.1, Kokkos 4.6.2, `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On
-D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off`, default Athena_M1_VGD_NMAX, no -DM1_VGD_NCCL. Build in a compute job.
-> `/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_1bd0c2ce_nofma`; its md5 = BQ_MD5 for all scripts below.
Production binaries (reference arm P): d505d21f nofma md5 1440c05320a7bef88c4ad960019b9840 (He giant), 5b7e1bc5 nofma md5
6252e58e5192ce65ff1f1f627e1ce90c (E1). Optional cross-check (5 min, debug): 1bd0c2ce default keys = 21c4839d on your E1S fin51 smoke deck.

## 2. Seed (login node; one COPY per case, used by every arm)
```
export CR2_X=/resnick/groups/carnegie_poc/jingze/rol2u_1015/caltech_rol2u BQ_MD5=<md5 of 1bd0c2ce>
bash $CR2_X/seed_cr2.sh hg ; bash $CR2_X/seed_cr2.sh e1
```
Takes the newest VALID rst (rst_pick.py) of hegiant_1014/wide4/run/rst and e1s_1014/prod/rst (>= 5 min old), `cp -p` into
rol2u_1015/<case>/seed/rst, cmp + md5, prints t / cycle / `<par_end>` offset. (If hgF4w link 4513229 is still pending, the He giant seed is
its final rst at t ~2.5918e6 = 30 d.) Edit only cr2_common.sh if a path/launch differs (LAUNCH_HG must equal the hgF4w launch: 1 rank per
H200, 4 ranks x 1 MeshBlock 897x64x64; LAUNCH_E1 = e1s_1014/launch2.sh, 2 ranks x 2 MeshBlocks 536x64x64). UCX_TLS=^xpmem is exported.

## 3. Gate (debug QoS, one job per case, <= 30 min; must PASS before any ROL2U long link runs)
```
sbatch -J cr2g_hg -A carnegie_poc -p gpu --qos=debug -t 00:30:00 --nodes=1 --ntasks=4 -c 8 --gres=gpu:nvidia_h200:4 --mem=240G \
       --exclude=hpc-sm-01-09 --export=ALL,CASE=hg,BQ_MD5=$BQ_MD5,CR2_X=$CR2_X $CR2_X/gate_cr2.sh
sbatch -J cr2g_e1 ... --ntasks=2 --gres=gpu:nvidia_h200:2 --export=ALL,CASE=e1,... $CR2_X/gate_cr2.sh   (after the hg one: debug MaxJobsPU 1)
```
Arms from the seed (hg 30 cycles, e1 20 cycles; hst every cycle, final rst/bins at the end): P production binary + production keys;
B = 1bd0c2ce + production keys + time/dt_retry=true (BASE); C = B + ROL2U; L1 (half) -> L2 (`-r` L1) ROL2U; F = C + forced retry at
seed cycle + 5 (fac 1.0). Gate lines (all must hold -> rol2u_1015/<case>/gate/VERIFIED):
- **GATE-BASE P = B BITWISE** (rst data after `<par_end>`, final bins via cmplast.py, all hst rows) = FIN4/5b7e1bc5 rst readable by
  1bd0c2ce and bitwise to production.
- **GATE-RESTART C = L2 BITWISE**; **GATE-RETRY C = F BITWISE** with exactly one `DT_RETRY ... reason=forced` line.
- rc 0, FATAL 0, NaN 0 in all six; `<par_end>` offset of the C and L2 rst < 44 000.
Report per arm: wall, s/cycle, Picard mean/max, NON-CONVERGED, lagged solves, floor clips, GPU peak (nvidia-smi).
If GATE-BASE fails: STOP, report the first differing hst row and the bin max rel diffs (do not run section 4/5). If only C = L2 or C = F fails:
report, hold the chains, run section 4 anyway (it is a test).

## 4. Short comparison window (normal QoS; after the gate; one job per case, BASE then ROL2U-window sequentially)
```
sbatch -J cr2w_hg -A carnegie_poc -p gpu --qos=normal -t 05:00:00 --nodes=1 --ntasks=4 -c 8 --gres=gpu:nvidia_h200:4 --mem=240G \
       --exclude=hpc-sm-01-09 --export=ALL,CASE=hg,BQ_MD5=$BQ_MD5,CR2_X=$CR2_X $CR2_X/win_cr2.sh
sbatch -J cr2w_e1 ... -t 06:00:00 --ntasks=2 --gres=gpu:nvidia_h200:2 --export=ALL,CASE=e1,... $CR2_X/win_cr2.sh
```
Window = seed t + 2.5e5 s (hg, ~13k cycles at dt ~19 s; B ~1.8 h) / seed t + 2e6 s (e1; wall-capped: B -t 3:30, RW -t 2:00 since the
post-ramp E1 runs ~3 s/cycle with retries; the analysis uses the common part). Production output cadence (rst_pick last_time keys), so
B's hst rows can be checked verbatim against the production run's own hst (HSTCMP below). RW = ROL2U over the same window in its own dir
(its end-of-run report gives BE steps / floor clips / lagged solves for exactly the window).
OPTIONAL (if GPU budget allows; own job): `ARMS=H` = BASE with cfl_number/2 (dt-error yardstick; hg -t 6:00, e1 -t 8:00, 4/2 H200).

## 5. ROL2U production-capable chains (submit HELD right after the gate job is queued; RELEASE only when <case>/gate/VERIFIED exists)
Same scripts, hooks, cadence and walltime policy as the production chains; only BIN and the ROL2U keys differ. tlim 1e9 (long-tlim rule),
QoS normal 7-d links (long QoS is taken by agcE1s, MaxJobsPU 1), 1 link + 2 afterany, athena -t = walltime - margin.
- **He giant `hgR2w`**, dir `rol2u_1015/hg/chain` (1 node x 4 H200, 4 x 897x64x64): copy hegiant_1014/wide4/run/run.cfg and
  link_hg4w_ext.sh; in run.cfg set BIN=.../athena_gpu_he_star_m1_1bd0c2ce_nofma, BIN_MD5=$BQ_MD5, TLIM=1e9,
  XKEYS="time/tlim=1e9 time/dt_retry=true <ROL2U keys>", UCX_TLS=^xpmem, and REMOVE any first-job seeding that picks/links a production
  rst (the wide4 run.cfg hard-linked the newest prod rst at job start): instead `mkdir -p chain/rst && ln seed/rst/<rst> chain/rst/`
  (a link of OUR copy only). STOP rules unchanged (rc != 0 / FATAL / NaN / DIVERGED (non-DT_RETRY) / NON-CONVERGED > 0), prune keep 3.
  rmax guard: a second rmax_guard.sh instance for this dir and job name hgR2w (RLIM 150, same as hgF4w).
  `sbatch --hold -J hgR2w -A carnegie_poc -p gpu --qos=normal -t 7-00:00:00 --nodes=1 --ntasks=4 -c 8 --gres=gpu:nvidia_h200:4 --mem=240G
  --exclude=hpc-sm-01-09 link_hg4w_ext.sh <chain dir>` + 2 x `--dependency=afterany:<prev>` (also held).
- **E1 `agcE1r`**, dir `rol2u_1015/e1/chain` (1 node x 2 H200, 4 blocks 536x64x64): copy e1s_1014/prod/run.cfg; set
  BIN=1bd0c2ce, BIN_MD5=$BQ_MD5, IN = the same deck 680dadb0 (tlim 1e9; dt_retry=true inside), KEYS="<ROL2U keys>", TLIM=1e9, STOP_NC=1,
  E1C_GATE=1, E1C_VERIFIED=rol2u_1015/e1/gate/VERIFIED, E1C_DIR=e1s_1014/scripts, LAUNCH=e1s_1014/launch2.sh, E1C_MARGIN=1200;
  seed: `mkdir -p chain/rst && ln seed/rst/<rst> chain/rst/`; after GATE PASS write chain/GO_e1 (text: "USER 10-15 via coordinator:
  ROL2U production-capable chain, gate <job> PASS"). `sbatch --hold -J agcE1r -A carnegie_poc -p gpu --qos=normal -t 7-00:00:00
  --ntasks=2 -c 8 --gres=gpu:nvidia_h200:2 --mem=240G --exclude=hpc-sm-01-09 link_e1.sh <chain dir>` + 2 afterany.
- Never more: no change to hgF4w / agcE1s; if a ROL2U chain STOPs, report and hold (no automatic restart).

## 6. Analysis (bundle: cr2ana.py = Raven drf.py adapted; bin_convert.py; cmprst.py; cmphst.py)
Window table (after section 4 [+ H]; R = the RW window arm, the chain is used only for long-term tables):
```
python3 -I $CR2_X/cr2ana.py --mode e1 --bc $CR2_X  rol2u_1015/e1/win/B R:rol2u_1015/e1/win/RW [H:rol2u_1015/e1/win/H]
python3 -I $CR2_X/cr2ana.py --mode hg --rph <R_ph at the seed time, cm, from the hgF4w guard/page> --bc $CR2_X \
        rol2u_1015/hg/win/B R:rol2u_1015/hg/win/RW [H:rol2u_1015/hg/win/H]
python3 -I $CR2_X/cmphst.py <prod dir>/<prefix>.user.hst rol2u_1015/<case>/win/B/<prefix>.user.hst <seed t>   (and .hydro.hst)
```
Prints: s/cycle (window cycles), Picard mean, DT_RETRY, NON-CONVERGED (lines + report), BE steps, E-floor clips (solved E <= e_floor;
written-back eint <= 0), lagged solves, FATAL; hst L_top/L0, L_bot/L0, KE_int, Picard window means + rms(arm - BASE) and slope; shell
metrics at common bin dumps: e1 rho/v_r at 1.05 R_ph (lrho105, vr105), vr110, rsurf, KEe; hg r12, rho/v_r at r12 and at R_ph, vmax.
`FLAG` lines = ROL2U - BASE rms above the cfl/2 - BASE rms (only with H). Expected from Raven AG Car drift: Ltop differences are at the
chaos floor in all arms (FLAG on Ltop alone with ratio ~1 is not a failure); Lbot / KE_int / lrho105 / KEe R-BASE well below H-BASE.
HSTCMP: BITWISE over the rows where the production run has reached (if production has not reached the window yet, say so).

Long-term (ROL2U chain vs the running FIN4 production, read-only), side-by-side at MATCHING model times, same tables as your e1s / e1w
and He giant page reports: E1 every ~1e6 s (pulse_ana.py: vrV vrM vr105 vr110 d105 d110 KEe rsurf L_top/L0 KE_int w_mlt Picard,
DT_RETRY, s/cycle); He giant every ~2 d model time (R_ph, r12, vmax, L_top, Mdot_top, s/cycle, DT_RETRY). Flag at once: FATAL/STOP,
DT_RETRY bursts, Picard mean > 20, infall < -5 km/s (E1), rsurf > 1.15 (E1), r12/R_ph trend departing from production by > 5 %.

## 7. Report: docs/handover/NOTE-2026-10-15-caltech-rol2u.md
md5 + HEAD of 1bd0c2ce; seed rst (path, t, cycle, md5, `<par_end>` offset) per case; gate table (6 arms, all lines, header offsets,
peaks); window table (cr2ana.py output + HSTCMP); chain job ids / held / released / first s/cycle; then the long-term side-by-side
tables as they come. Any non-bitwise gate line in full (first differing hst row, max rel diff per variable).

## GPU budget (H200-h, approx.)
gate: hg 2 + e1 1 (debug) | window B+RW: hg ~19, e1 ~11 | optional H: hg ~22, e1 ~14 | chains: hg 672 per 7-d link, e1 336 per link.
