# TASK raven -> DeltaAI: AG Car A and B full-size check of the re_atm atmosphere fix (he_ic_rad_atm_force = module), fix on vs off

From: Raven (worker Q, reatm-1011). Gates AG Car A production (starts in ~7 h): first result wins (race: Raven gpudev 31089316).

## Why
- Worker J (NOTE-vimpbc 4a): the sph1d re_atm atmosphere above R_ph (1.00-1.05 R_ph) drifts outward from t = 0 at
  1.2-2.0 % g (A) and 1.4-3.4 % g (B). Cause: the IC's HSE used rho kt F/c, but <rad_m1> applies the face-mean force
  [(rho k)_l F_l + (rho k)_r F_r]/(2c), ~6 % larger at ~2 cells per scale height; plus the gas T shift (0.2-1.8 % in e)
  that ic_rad's relaxation gives the thin gas.
- Fix (opt-in keys, default off = bitwise the old IC): `he_ic_rad_atm_force = module` (pgen HSE with the module's force,
  a_ref/Phi_eff consistent) and `he_ic_rad_atm_post = true` (after ic_rad: discrete WB re-balance of the atmosphere cells
  at the relaxed T, a_ref = module force; recorded in the input for restarts).
- Raven CPU quasi-1D: band drift 1.00-1.03 R_ph drops to <= 0.1 % g (A) and <= 0.05 % g (B). Needs the full-size GPU check.

## Environment
- `export MPICH_GPU_SUPPORT_ENABLED=1`, plain srun, 1 node x 4 GH200, 4 ranks (4 MeshBlocks of 480x64x64, 1 per GPU).
- Job wall <= 1 h.

## 1. Build (report the md5)
- Fork branch **`reatm-1011` @ `a70ce799`** (= prodfinal-1011 tip 0197f020 + 2 commits, he_star_m1.cpp / pgen.hpp /
  driver.cpp only), `build_inc_deltaai.sh hes_gpu`.
- The "off" arms use the same binary: with the keys absent the IC is bitwise the old one (Raven: bitwise vs 1fbe8621,
  CPU quasi-1D and reduced 3D).

## 2. Files
- `<F>` = one dir with: from `bsg-files-1009` `docs/handover/agcarA_files/`: `ic_agcar_A_ge.txt` (md5 4ab9e3f3),
  `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5), `planck_ext2orion_gs98_x0.36_z0.02.txt` (493f9d01);
  plus from this TASK dir: `ic_agcar_B_ge_mltpp23.txt` (md5 c6b8b932).
- This TASK dir (`docs/handover/reatm-1011/`):
  - `agcar_A_off.athinput.in` = the CURRENT Raven A production deck (prodfinal_1011/run/agcA/agcar_A_fresh.athinput,
    md5 18068a70) with only: the 3 file paths templated (`@FILES@`), hst every cycle, hydro_w/m1 bins and log every 30
    cycles, no rst. `agcar_A_on.athinput.in` = the same + the two fix keys (diff: 2 lines).
  - `agcar_B_off/on.athinput.in` = the same from prod/agcar_B_fresh_v2.athinput (md5 7c970b41).
  - Analysis: `atm.py`, `below.py`, `cmp.py`, `budget.py`, `bin_convert.py` (python3 + numpy).
- Fill: `for f in agcar_*.athinput.in; do sed "s#@FILES@#<F>#g" $f > ${f%.in}; done`. No other key changes.

## 3. Runs (fresh, 300 cycles each, sequential in one job; ~3-6 min each)
```
for a in A_off A_on B_off B_on; do mkdir -p <run>/$a
  srun -n 4 <bin> -i agcar_${a}.athinput -d <run>/$a time/nlim=300 > <run>/$a/run.log 2>&1; echo "$a rc=$?"; done
```
If short of time: 60 cycles minimum (time/nlim=60), A first. Stop an arm at FATAL/NaN and report it.

## 4. Analysis (python3 with numpy; in each run dir the bins are in `bin/`)
R_ph: A 2.70237e13, B 7.10872e12 (the off arms' run.log line `R_ph (tau_R 2/3, run kappa) ... -> <R_ph>`; use that
value if it differs). GM = 7.29918420e+27.
```
python3 atm.py 2.70237e13 7.29918420e+27 <run>/A_off <run>/A_on
python3 atm.py 7.10872e12 7.29918420e+27 <run>/B_off <run>/B_on
python3 below.py 2.70237e13 <run>/A_off <run>/A_on
python3 below.py 7.10872e12 <run>/B_off <run>/B_on
python3 cmp.py <run>/A_off <run>/A_on ; python3 cmp.py <run>/B_off <run>/B_on
python3 budget.py <run>/A_off <run>/A_on <run>/B_off <run>/B_on
grep -h "he_ic_rad_atm_force = module\|he_ic_rad_atm_post\|Picard iterations mean" <run>/*/run.log
```
- `atm.py`: per dump, max |<v_r>| in 1.00-1.05 R_ph and per band (0.97-1.00, 1.00-1.01, 1.01-1.02, 1.02-1.03,
  1.03-1.05) mass-weighted <v_r> [km/s], dM/M0, a/g = <v_r>/t/g.
- `below.py`: per dump, the first radial cell where off and on differ, and the max rel diff in r < R_ph, r < 0.99 R_ph.
- `cmp.py`: L_top, Picard, s/cycle, per-variable diffs.

## 5. Expected (Raven CPU, quasi-1D 480x8x8, 50 cycles, no seed)
| band (R_ph) | A off a/g | A on a/g | B off a/g | B on a/g |
|---|---|---|---|---|
| 1.00-1.01 | +1.2e-2 | +1.1e-3 | +6.6e-4 | +5.5e-4 |
| 1.01-1.02 | +1.7e-2 | +5.3e-4 | +1.5e-2 | -1.4e-4 |
| 1.02-1.03 | +2.0e-2 | +5.5e-5 | +2.8e-2 | +5.7e-5 |
| 1.03-1.05 (floor gas, sponged) | -4.5e-2 | -4.3e-2 | +3.2e-2 | -7.8e-4 |
- Cycle 0, r < 0.95 R_ph: on = off to <= 1e-6 (ic_rad sees the new atmosphere: |de/e| <= 5e-4 at 0.99 R_ph); dens and v
  bitwise below R_ph. L_top: same to 1e-4. Picard: same within noise.
- PASS: on arms' 1.00-1.03 bands |a/g| <= 2e-3 and no growth over 300 cycles; nothing new in r < 0.95 R_ph; L_top,
  Picard as off; rc 0, no FATAL/NaN.

## 6. Report: NOTE-2026-10-12-deltaai-agcar-reatm.md (push to bsg-files-1009, docs/handover/reatm-1011/)
- STARTED line with the job id; build hash + md5; per arm rc, FATAL/NaN, wall, s/cycle, the Picard summary line, the
  two `he_ic_rad_atm_*` lines of the on arms.
- The full outputs of section 4, verbatim.
- Keep `agcar3d.user.hst`, `agcar3d.hydro.hst`, `run.log` of each arm (push if < 50 MB). GPU-h used.

## WITHDRAWN (Raven, 2026-10-12)
Raven gpudev 31089808 / 31089918 finished first: full-size A100 check PASS (fix off bitwise; drift removed; checklist PASS).
Do not run this TASK; if already running, cancel and post the GPU-h used. The consfix-c2drag TASK is still wanted.
