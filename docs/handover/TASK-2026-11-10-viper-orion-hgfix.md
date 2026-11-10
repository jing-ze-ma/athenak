# TASK 2026-11-10 viper -> ORION: hgfix -- CPU gates + small-state repro of the HG stage-2 eint < 0 lock

User GO 13:50 PST ("if GPUs here are not enough, use Orion for tests"). CPU work only. Answer as NOTE hgfix on
bsg-files-1009 ("NOTE back"). The full-HG A/B stays on GPU (viper 6-node successor / Raven), not here.

## Why
The HG (He giant, cubed sphere) production runs in permanent backward Euler: 0 of 812 hesdirk2 stage steps succeed
(viper hgdiag_1110b AB2). Every retry fails at stage 2 in near-vacuum fast cells (r/R_ph 0.82-0.90, rho/<rho> ~1e-9,
|v| ~1500 km/s): the gas old vector of the stage solve has eint = eint_h + t2inc_EN - dKE < 0 because the momentum
increment t2inc_m (from the stored slope K_m) carries kinetic energy that K_EN does not (gas split: cK_KE +0.145
cancelled by post -0.145, dKE 0.074 = 1e6 eint_h). Candidates for the momentum without energy, both after the
Picard solve and both paying their KE from e_int: the implicit_am_cs L deposit (rad_m1_implicit.cpp m1_am_cs) and
the conservative write-back's F-clip momentum (rad_m1_impl_wb_kernel.hpp, consrc-1011 clip block).

## Code
Branch **fork/hgfix-1110** (fork/hgkmom-1110 a59e1032c + fork/prod-1113-cand 33312dd7d merged + the fix). Commit
**3d4fd7da5**. New opt-in keys (read only when named; absent = bitwise):
- `<rad_m1>/time2_inc_kecap = th` (0..1): at both stage formations, where eint(U_h + t2inc) < th eint(U_h) the cell's
  t2inc momentum is scaled down (total energy kept, dropped momentum counted; INCKECAP lines).
- `<rad_m1>/implicit_clip_kecap = f` (0..1): the F-clip momentum gives the gas at most f e_int of KE (state + slope).
- `<rad_m1>/implicit_am_cs_kecap = f` (from hgkmom-1110): the am_cs deposit gives at most f e_int (state + slope).
- DEBUG `<problem>/dbg_vac_seed = f` + `dbg_vac_rlo`, `dbg_vac_rhi` (code units = cm), `dbg_vac_frac` (0.05),
  `dbg_vac_vr` (cm/s, 0), `dbg_vac_omega` (rad/s, 0; cubed sphere): fresh start only; every cell of the shell first
  gets a solid-body rotation about z, then the cells whose hash < frac get gas x f (T, v kept) and v_r.

Build (as NOTE-2026-11-08-orion-viper-consfix-ledger):
```
git -C /orion/u/jinma/ATHENAK/athenak fetch origin
git -C /orion/u/jinma/ATHENAK/athenak worktree add --detach /orion/ptmp/jinma/Athenak/wt_hgfix1110 <commit>
cd /orion/ptmp/jinma/Athenak/wt_hgfix1110 && git submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build_hf -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_hf -j32
```
Reference binary (keys-off gate): the same recipe at **fork/prod-1113-cand 33312dd7d** into build_ref.

## Inputs (on this branch)
- `docs/handover/hgfix-1110/hg_small_deck.athinput`: the HG production deck on a small cubed sphere (704 x 16 x 16
  blocks, 6 MeshBlocks, implicit_am_cs + am_transport, hesdirk2).
- `docs/handover/hgfix-1110/prodkeys.txt`: the production command-line keys (one line), appended to every run.
- `docs/handover/hgfix-1110/mkdeck.py IN OUT block/key=v ...`: inserts keys into their block (a fresh start refuses
  command-line keys the deck lacks).
- Tables: `docs/handover/hgpos-1109/files/` (ic_hegiant_14p5M_mk.txt, rosseland_mesakap_afe2dd8da2bd.txt,
  planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt; MD5SUMS there).

Common args (T = the tables dir):
```
B="problem/he_ic_file=$T/ic_hegiant_14p5M_mk.txt problem/he_opac_table=$T/rosseland_mesakap_afe2dd8da2bd.txt \
problem/he_planck_table=$T/planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt output1/dt=1e-30 output2/dt=1e30 \
output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 mesh/nx2=16 mesh/nx3=16 $(cat prodkeys.txt)"
srun -n NP --exact --mem=15G --cpu-bind=cores BIN -i deck.athinput $B time/nlim=NL -d RUNDIR
```
(deck.athinput = mkdeck.py hg_small_deck.athinput deck.athinput <arm keys>; arm keys go into the deck, not on the
command line.)

## Runs and pass criteria
A. Bitwise, keys absent, NL 4: ref vs new at NP 6, 3, 1. Pass: hegiant.hydro.hst and hegiant.user.hst identical
   (cmp) for each NP. Also new NP 6 vs new NP 3 with ALL keys on (below): identical.
B. Keys on, no seed, NL 4, NP 6: ALL = `rad_m1/implicit_am_cs_kecap=0.5 rad_m1/implicit_clip_kecap=0.5
   rad_m1/time2_inc_kecap=0.5`. Pass: no FATAL; report whether hst is identical to ref (it is if no cap binds:
   report the kecap / INCKECAP counters).
C. Repro, NP 6, NL 120, SEED = `problem/dbg_vac_seed=1e-9 problem/dbg_vac_rlo=6.52e12 problem/dbg_vac_rhi=7.15e12
   problem/dbg_vac_frac=0.05 problem/dbg_vac_vr=1.5e8 problem/dbg_vac_omega=2e-7`;
   arms S0 / SA / SC (= SEED + implicit_clip_kecap 0.5) / SI / S3 as below. arms S0 = SEED, S1 = SEED + am_cs_kecap 0.5,
   S2 = SEED + time2_inc_kecap 0.5, S3 = SEED + ALL. Report per arm: rc, FATAL count, `NOT ADMISSIBLE` count (stage 1/2),
   the end summary lines `time2_auto_be=...skipped-stage BE steps ... stage retries (succeeded)`, `wall split`, the
   kecap / clip_kecap / INCKECAP counters, and hydro.hst totals (mass, energy) at the last common row vs S0.
   Pass for the repro: S0 shows stage-2 NOT ADMISSIBLE / auto_be steps (the lock), and S3 shows fewer of them with
   its caps binding (counters > 0); conservation: mass identical, total energy change explained by the counters.
   Viper measured (login, throttled): seed without rotation, 12 cycles: 11 stage steps, 0 faults; seed + omega 2e-7,
   10 cycles: 9 stage steps, 0 faults, no cap bound (am_cs deposits sum |L_z| 9.5e43 vs 8.3e40 without rotation), keys
   on bitwise with keys off.  So the early cycles do not lock: run NL 120 (and 400 if S0 is still clean), and if S0
   stays clean add a second pair with dbg_vac_seed=1e-11, dbg_vac_frac=0.2.  The same repro runs on GPU (Raven, mesh
   N 32, 200 cycles); do not duplicate it on GPU.

## NOTE back
NOTE-2026-11-1x-orion-hgfix.md on bsg-files-1009 with: commit built, md5 of both binaries, the tables of A/B/C, the
counter lines verbatim (cut to 300 chars), and anything that FATALed with its last 20 log lines.

Commit: fork/hgfix-1110 3d4fd7da5 (reference fork/prod-1113-cand 33312dd7d).
