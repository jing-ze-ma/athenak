# TASK raven -> DeltaAI: full-size BSG timing of the bsgaccel-1011 speedups (1 short job, ghx4-interactive)

From Raven worker M (bsgaccel-1011). Raven cannot run the full 256^3 BSG on one gpudev node (A100 40 GB). This job
confirms the full-size speedup on GH200 that Raven measured on the reduced BSG (1.8-1.97x, bsgaccel NOTE docs/handover/NOTE-2026-10-12-raven-bsgaccel.md). Limits: ≤ 40 min wall,
1 node x 4 GH200, 1 rank per GPU, plain srun, `export MPICH_GPU_SUPPORT_ENABLED=1`. Reply SKIP if the queue is
more than about 2 h.

## 1. Builds (`build_inc_deltaai.sh hes_gpu`, PROBLEM=he_star_m1; report the md5s)
- BASE: fork `prodfinal-1011` @ **422c2af7** (the tip before the bsgaccel merge).
- NEW: fork `prodfinal-1011` @ **0197f020** (= 422c2af7 + merge of bsgaccel-1011: 3 bitwise speed commits).
  If the tip has moved past 0197f020 with code changes, still build 0197f020.

## 2. Input
`bsgR` from your bsgrec smoke (docs/handover/bsgrec-1011/bsgR.athinput.in, files dir resolved). Make one copy with
these lines added in `<rad_m1>`, all named at their defaults so that the command-line keys exist:
```
implicit_precond_float = false
vet_gd_halo_ipc = false
implicit_mg_levels = 3
implicit_mg_halo = true
implicit_gf_modes = 2
implicit_cs2_shift = false
```
`implicit_precond` and `vet_gd_twin_lowmem` are already in the deck.

Common keys: `time/nlim=24 time/ndiag=1 output1/dt=1e-6 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30
output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`.

## 3. Runs (in this order; stop on FATAL)
| run | binary | ranks x blocks/GPU | extra keys |
|---|---|---|---|
| B4 | BASE | 4 x 4 | - (deck: lowmem true) |
| N4 | NEW | 4 x 4 | `rad_m1/implicit_precond=mg_gf rad_m1/vet_gd_halo_ipc=true` |
| N4n | NEW | 4 x 4 | N4 + `rad_m1/vet_gd_twin_lowmem=false` (OOM is an answer) |
| B2 | BASE | 2 x 8 | - |
| N2 | NEW | 2 x 8 | as N4 |
| C4 | NEW | 4 x 4 | `rad_m1/implicit_gas_src=conservative2` (lowmem true, default precond) |
| C4g | NEW | 4 x 4 | C4 + `rad_m1/implicit_precond=mg_gf rad_m1/vet_gd_halo_ipc=true` |

**C4 / C4g are URGENT (run them FIRST if time is short):** the old `conservative` gas source blew up on this BSG in your
bsgrec smoke (RC4: stage 1 NOT ADMISSIBLE, dt 87 -> 0.008 s). Production now uses its fix `conservative2`; we need to know
whether it is healthy on the full BSG.

## 4. Report: NOTE-2026-10-12-deltaai-bsgaccel-fullbsg.md
- Per run: rc, FATAL, NON-CONVERGED, the Picard line, the bicgstab `inner iterations mean` line, peak MiB per GPU.
- s/cycle from the `elapsed=` lines: steady (cycles 3..22 that are not multiples of 10), the rotation cycles 10 and
  20, and the average.
- The line `vgdrot-1011 vet_gd_halo_ipc: ON ...` must be present in N runs.
- B4 vs N4 at cycle 24, from user.hst: max relative difference of L_int, L_top, E_rad, Etot over the rows.
  Expected ≤ 1e-7.
- C4 / C4g: dt per cycle (from user.hst or the cycle lines), the `time_scheme=hesdirk2: stage steps= ... backward-Euler steps= ...
  stage fallbacks=` line, any `NOT ADMISSIBLE` lines, vimp positivity fallbacks, min E, gas Newton fallbacks, L_top per cycle
  (paste the hst rows), compared with your R4 (matched) and RC4 (old conservative) from the bsgrec smoke.
- GPU-h used.
