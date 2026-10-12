# TASK raven -> DeltaAI: GPU cost of vet_gd_rotate_slope = correct on the full AG Car A deck (keep vs correct)

From: Raven (worker R2, rotslope-1012). This is a cost check only, not production. Raven gpudev has about 10 jobs queued
(estimated wait over 2 h). Raven job 31105120 races this one: the first result wins.

## Why
- With the physical-time vet_gd rotation (vet_gd_rotate_dt), the hesdirk2 FSAL slope K1 of the first step of a new
  rotation block belongs to the OLD quadrature set. That leaves a first-order time error in thin cells.
- New opt-in key `vet_gd_rotate_slope = correct`: on each such step there is one extra stage-type implicit solve
  (correction solve) with the new set. The vet_gd sweep is reused, so there is no extra sweep.
- Expected overhead: one Picard solve per rotation. AG Car A rotates every 20 cycles (vet_gd_rotate_dt 2.5e4, dt 1250).
- Question: the s/cycle of correct vs keep (target: at most a few %), and that rc, Picard and NON-CONVERGED are fine.

## Environment
- `export MPICH_GPU_SUPPORT_ENABLED=1`, plain srun, 1 node x 4 GH200, 4 ranks (4 MeshBlocks of 480x64x64, 1 per GPU).
- Job wall <= 1 h (expected ~15 min).

## 1. Build (report the md5)
- Fork branch **`rotslope-1012` @ `1bfe213a`** (= method-1012 e39ce690 + 1 commit: driver.cpp, rad_m1.cpp/.hpp,
  rad_m1_implicit.cpp, rad_m1_time2.cpp, rad_m1_vetgd.cpp). Build with `build_inc_deltaai.sh hes_gpu` (he_star_m1 pgen).
- Both arms use the same binary. They differ only in the key `vet_gd_rotate_slope` (keep / correct).

## 2. Files
- `<F>` = one dir with, from `bsg-files-1009` `docs/handover/agcarA_files/`:
  - `ic_agcar_A_ge.txt` (md5 4ab9e3f3)
  - `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5)
  - `planck_ext2orion_gs98_x0.36_z0.02.txt` (493f9d01)
- This TASK dir (`docs/handover/rotslope-1012/`):
  - `agcar_A_FIN.athinput.in` = the method-1012 FIN full AG Car A deck (Raven final_1012/gpu/inp/agcA_full_FIN.athinput,
    480x128x128), with only the 3 file paths templated (`@FILES@`).
  - `agcar_A_FINC.athinput.in` = the same deck with `vet_gd_rotate_slope = correct` (a 1-line diff).
  - `tana.py`: per-cycle timing from the `elapsed= cycle=` lines.
- Fill the templates (no other key changes):
  `for f in agcar_A_*.athinput.in; do sed "s#@FILES@#<F>#g" $f > ${f%.in}; done`

## 3. Runs (fresh, 120 cycles each, sequential in one job, FINC first)
```
O="time/ndiag=1 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output5/dt=1e30 time/nlim=120"
for a in FINC FIN; do mkdir -p <run>/$a
  srun -n 4 <bin> -i agcar_A_${a}.athinput -d <run>/$a $O > <run>/$a/run.log 2>&1; echo "$a rc=$?"; done
```
- If short of time, run 60 cycles minimum.
- Stop an arm at FATAL or NaN and report it.

## 4. Analysis
```
for a in FINC FIN; do echo "$a $(python3 tana.py <run>/$a/run.log 2 20)"; done
grep -h "Picard iterations mean\|time_scheme=hesdirk2:" <run>/*/run.log
grep -h "cycle=" <run>/FINC/run.log | tail -1; grep -h "cycle=" <run>/FIN/run.log | tail -1
```
- `tana.py run.log 2 20`: steady = cycles with c % 20 != 0, rot = the rotation cycles (c % 20 == 0), avg = all
  cycles from 2 on.
- In the FINC log, the `time_scheme=hesdirk2:` line ends with `vet_gd_rotate_slope=correct solves=N failed=M`.
  Expected N = 6 for 120 cycles, M = 0.

## 5. Expected (Raven CPU reduced 480x32x32, 300 cycles)
- FINC: rc 0, Picard mean 3.23, max 6, 0 NON-CONVERGED, 14 correction solves, 0 failed. FIN: Picard 3.15.
- Raven A100, FIN on 78b2c961: 0.335 s/cycle average, 0.331 steady, 0.45 on rotation cycles.
- Expected for FINC: rotation cycles about +1 Picard solve, i.e. roughly +0.1 to 0.2 s on A100. The average should be
  at most a few % above FIN.
- PASS: both rc 0, no FATAL or NaN, 0 NON-CONVERGED; FINC avg / FIN avg <= 1.05.

## 6. Report: NOTE-2026-10-12-deltaai-rotslope-cost.md (push to bsg-files-1009, docs/handover/rotslope-1012/)
- STARTED line with the job id; build hash and md5.
- Per arm: rc, FATAL/NaN, wall time, the tana.py line, the Picard line, the hesdirk2 line.
- Keep `run.log` of each arm (push it). GPU-h used.
