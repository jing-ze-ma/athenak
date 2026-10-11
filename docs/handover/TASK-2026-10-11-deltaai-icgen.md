# TASK raven -> DeltaAI: AG Car A, transport-consistent radiation IC `he_ic_rad = sph1d` (1 short ghx4-interactive job)

From the Raven session (icgen-1011). This is a short test, matching DeltaAI's role (NOTE-2026-10-09-deltaai-role-short-tests).
The same arms are queued on Raven gpudev as job 31068342. **First result wins**: post your NOTE as soon as you have it.
If ghx4-interactive will not start within ~2 h, reply SKIP in the NOTE.

## Background (already measured on Raven; do NOT repeat)
`problem/he_ic_rad = sph1d` solves the 1-D spherical grey RT problem in radiative equilibrium (VET: p-z formal solution
plus moment equation) on the IC column at setup, with the run's own kappa_t. It replaces E above tau_R = 3000.
Optionally it also sets gas T (radiative equilibrium) and re-integrates HSE above 0.9 R_ph.

Raven gpudev job 31067719 (60 cycles, full size), gas = re + hse:
- L_top step-1 spike +53.6 % -> +1.0 %.
- E/J at the first build: 1.06 -> 0.018.
- The photosphere stays static (0.97-1.028 R_ph: max |v_r| 0.14 km/s, rho drift 0.9 % at cycle 50; default 2.35 km/s, 33 %).
- Runs are bitwise run to run.

The open issue is cost. With gas = re, the implicit solve takes ~6 Picard iterations per step from cycle 60 on, against
~3.2-4 for the default. The cause is the cold floor gas (rho = 1e-16 above 1.03 R_ph) falling 3x faster. The new variant
`he_ic_rad_gas = re_atm` (radiative equilibrium only above the floor) is meant to fix that. **This task tests re_atm.**

## Code (fork github.com/jing-ze-ma/athenak), PROBLEM=he_star_m1, `build_inc_deltaai.sh hes_gpu`
- Branch `icgen-1011`, commit **3d8861c7** (on icrad-1011 7045bc14; the default path is bitwise the old IC: checked
  against vetacc j.31064135 over 60 cycles on A100).
- Job script: `export MPICH_GPU_SUPPORT_ENABLED=1` (without it, Cray MPICH fails at startup with process_vm_readv).

## Files
`docs/handover/agcarA_files/` on bsg-files-1009: the IC, the Rosseland and Planck tables, and the deck agcar_A_va.athinput.
Point he_ic_file / he_opac_table / he_planck_table at your copies; md5s as in TASK-2026-10-11-caltech-agcarA-face-icrad.md.
No new data file is needed (the solve runs inside the pgen).

## Runs: one job, 1 node x 4 GH200, 1 rank per GPU (4 blocks 480x64x64), --time 1:00, plain srun
Common keys in every run (set them in a per-run copy of the deck, or on the command line if your deck names them):
`time/nlim=120 time/ndiag=1 rad_m1/implicit_cfl=1251.3` (dt = 1239.996 s), hst every cycle (deck default),
`output2/dcycle=10 output3/dcycle=10 output7/dcycle=10 output6/dcycle=100000`, no restart output.
Also add `rad_m1/ic_rad_diag=2` to every run (read-only; prints `ic_rad run_c0 ... max|<E>/<J>-1|`).

Key sets (add to `<problem>` / `<rad_m1>`):
- SPH  = `problem/he_ic_rad=sph1d problem/he_ic_rad_gas=re_atm problem/he_ic_rad_hse=true`
- SPHR = `problem/he_ic_rad=sph1d problem/he_ic_rad_gas=re problem/he_ic_rad_hse=true` (the old variant, for the Picard comparison)
- FR   = `rad_m1/ic_rad=fs_relax rad_m1/ic_rad_relax_n=40 rad_m1/ic_rad_relax_tol=1e-4`

| arm | extra keys |
|---|---|
| g0   | (none: the default IC) |
| gha  | SPH |
| gha2 | SPH (repeat: the run-to-run check) |
| ghafr | SPH + FR |
| gra  | `problem/he_ic_rad=sph1d problem/he_ic_rad_gas=re_atm` (no hse) |
| gh   | SPHR |

## Report: NOTE-2026-10-11-deltaai-icgen.md
- Every run: rc, FATAL, NON-CONVERGED, Picard mean/max (the end-of-run line), and s/cycle (median, excluding rotation
  cycles 10, 20, ... and the last). Also the setup line `he_ic_rad = sph1d (gas ...) ... N s` (the solve cost).
- **Picard per cycle**, as means over cycles 1-20, 21-40, ..., 101-120 (hst column Picard), for every arm. This is the
  main question: is gha at the g0 level (about 3.2-4) or at the gh level (about 6)?
- L_top per cycle 0-10, and at 20, 60 and 120, normalised to g0 cycle 0.
- The `ic_rad run_c0` line (max |<E>/<J>-1|, top value) per arm, and for ghafr the `ic_rad: done` line.
- Photosphere hydrostatics at cycles 10, 30, 60 and 120: max |v_r| and max |<rho>/<rho>(cycle 0) - 1| over 0.97-1.028 R_ph
  (R_ph = 388.3 Rsun = 2.7014e13 cm). Also the shell mean v_r at 2.0 R_ph (the floor-gas infall).
- gha vs gha2: bitwise? (all bin files and the hst)
- GPU-h used.
