# TASK raven -> DeltaAI: AG Car A follow-up, ic_rad restart safety + hesdirk2 + halo (1 short ghx4-interactive job)

Thanks for NOTE-2026-10-11-deltaai-agcarA-face-icrad: it answered the main questions first. This is the small rest.
Raven has it queued too (31066834 gpudev / 31066835 gpu); **first result wins**. Same binaries, files and environment as
that task (remember MPICH_GPU_SUPPORT_ENABLED=1). 1 node x 4 GH200, --time 0:30. Fixed dt implicit_cfl=1251.3 unless noted.

## Runs
Binary B (icrad-1011 7045bc14):
1. `berx`: ic_rad=relax (ic_rad_relax_n 40, ic_rad_relax_tol 1e-4), 20 cycles, with restart output every 10 cycles
   (an rst at cycle 0 just after the IC step if the code writes one, else at cycle 10), hst every cycle.
2. `rs10`: restart from berx's cycle-10 rst and run to cycle 20. **Check: the IC step must NOT run again on restart** (no
   relax/fs messages in the log), and cycles 11-20 must be bitwise identical to berx's (hst + final bin).
   If a cycle-0 rst exists, do the same from it (`rs0`, to cycle 10).
3. `hprx`: ic_rad=relax with `rad_m1/time_scheme=hesdirk2 rad_m1/time2_vet_col=predict rad_m1/implicit_vimp=true
   rad_m1/time2_vstage=true`, 60 cycles. Compare L_top(t) with your B relax run.
Binary A (deepatm-1011 e84d2513):
4. `hc1`: vet_bc_diff=face with `rad_m1/vet_gd_halo_compact=1` (production default is 2), 50 cycles. Is it bitwise vs your
   `A face` run (bin + hst)?

## Report: NOTE-2026-10-11-deltaai-agcarA-icrad-restart.md
rc / FATAL / NON-CONVERGED / Picard / s/cycle per run; restart: IC step re-run yes/no, and bitwise yes/no (which files);
hprx: L_top table cycles 0-10, 30, 60 next to your B relax; hc1: bitwise yes/no. GPU-h.

## Addendum (raven, user request): does relax perturb the envelope gas?
5. From your existing NOTE-2026-10-11-deltaai-agcarA-face-icrad runs (B default vs B relax / fs_relax, dt), no new runs:
   - at cycle 0 (after the IC step): the max relative change, relax vs default, of the shell-mean gas eint (or T) and
     the gas pressure, in the bands < 0.8, 0.95-1.05, 1.05-1.5 and > 1.5 R_ph;
   - over cycles 0-60: KE_int and KEr_int from the hst, and the max |velx| per dump, default vs relax / fs_relax.
   Say whether relax creates a velocity transient (and how large compared with the default IC's).
