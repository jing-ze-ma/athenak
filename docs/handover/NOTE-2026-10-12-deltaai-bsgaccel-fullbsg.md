# NOTE DeltaAI: full-size BSG timing of bsgaccel + conservative2 health check (TASK-2026-10-12-deltaai-bsgaccel-fullbsg) -- DONE

**Picked up** 2026-10-11 19:50 CDT. Building NEW 0197f020 then BASE 422c2af7. Deck = bsgR (bsgrec smoke) + the 6 named-default
`<rad_m1>` keys. Order in one job: **C4, C4g first**, then N4, N4n, N2 (NEW), then B4, B2 (BASE). Job id + STARTED line follow.

- 19:57 CDT: NEW built (athena_hes_gpu_0197f0202e10, md5 d36ed4ff9aed3bc092ab03e42d425c93); BASE 422c2af7 building (job waits up to 10 min for it after the NEW runs).
  Job **3363179** PENDING, priority 2173, Slurm start estimate **21:53 CDT (~2 h)**: at the TASK's SKIP limit; I keep it queued unless you post WITHDRAWN.

**STARTED**: job 3363179 RUNNING on gh152 since 20:09:04 CDT (C4 first). Raven: no need to run this.

## URGENT, C4 (conservative2, default precond, lowmem): NOT healthy on the full BSG (posted 20:13 CDT, job still running)
rc 0, FATAL 0, but **NON-CONVERGED = 6**, Picard mean 34.9 / max 200 (48 solves), and hesdirk2 stage 1 is **NOT ADMISSIBLE on every
cycle 2-23** (`stage steps=1 backward-Euler steps=23 stage fallbacks=22`); vimp positivity fallbacks 41, **min E from the solve
-4.05e6**; gas Newton fallbacks 4.5e8 (1.6e-2 per cell-pass; R4 matched 1.2e-3).
**dt collapses** at cycle 3 from 87.37 s to 0.00836 s, doubles back to 2.14 s by cycle 11, collapses again to 0.0015 s at cycle 12
and to 0.019 s at cycle 19 (t only reaches 266.8 s in 24 cycles vs 2097 s at full dt).
L_top: 1.6716e38 (c1) 4.457e37 (c2) 1.854e37 (c3), then 1.77-2.69e37 (R4 matched: 4.32-4.75e37 at cycles 3-10). Unlike the old
`conservative` (RC4: L_top -> 8.7e31) it does not lose the luminosity, but it does not take a single hesdirk2 step after cycle 1.
Identical to R4 / RC4 through cycle 1 (L_top 1.67165e38). Peak 49 789 MiB per GH200 (4 blocks/GPU). C4g and the timing runs follow.

## FINAL (20:20 CDT): job 3363179 complete, 7/7 runs rc 0, FATAL 0 (gh152, wall 4 min 13 s)
Binaries: NEW 0197f020 `athena_hes_gpu_0197f0202e10` md5 d36ed4ff9aed3bc092ab03e42d425c93;
BASE 422c2af7 `athena_hes_gpu_422c2af767be` md5 5d8f8a71c499992b12a9fa9dc58facc9.
Deck = bsgR + the 6 named-default keys; common command-line keys as in the TASK (all outputs off except hst every cycle).

### Summary
- **Full-size speedup NEW (mg_gf + halo_ipc) vs BASE, steady s/cycle: 1.60x at 4 x 4** (0.665 -> 0.415 mean; median 0.645 -> 0.409)
  and **1.79x at 2 x 8** (1.013 -> 0.565 mean; median 1.010 -> 0.539). Average incl. rotation and last cycle: 1.48x / 1.59x.
- **lowmem=false gives nothing** at 4 blocks/GPU (N4n 0.421 vs N4 0.415 s/cycle) and costs 79.2 GB vs 49.4 GB per GPU.
- **B4 vs N4 hst max rel diff: L_int 6.0e-8, L_top 3.9e-10, E_rad 3.5e-13, Etot 1.0e-13** (<= 1e-7 as expected; not bitwise because the
  preconditioner differs). N4 = N4n = N2 identical solver counters (Picard 2.72 / 7); B4 = B2.
- **conservative2 is NOT healthy on the full BSG, with either preconditioner** (C4 and C4g: NON-CONVERGED 6 / 4, stage 1 NOT ADMISSIBLE
  on every cycle 2-23, dt collapse to ~1e-3-1e-2 s). Matched (N4/B4) is clean on the same deck and binary.

### Per run
| run | bin | layout | NON-CONV | Picard mean / max | bicgstab inner mean / max | stage steps / BE / fallbacks | min E | s/cycle steady mean (median) | c10 / c20 | average | peak MiB/GPU |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C4 | NEW | 4x4 | **6** | 34.9 / 200 | 1.86 / 85 | 1 / 23 / 22 | **-4.05e6** | 1.782 (1.119) | 4.07 / 0.60 | 1.903 | 49 789 |
| C4g | NEW | 4x4 | **4** | 31.7 / 200 | 1.02 / 35 | 1 / 23 / 22 | **-2.97e4** | 1.801 (1.250) | 4.45 / 1.20 | 2.029 | 50 219 |
| N4 | NEW | 4x4 | 0 | 2.72 / 7 | 0.98 / 2 | 23 / 1 / 0 | 7.87 | **0.415** (0.409) | 0.350 / 0.351 | 0.522 | 49 391 |
| N4n | NEW | 4x4 lowmem off | 0 | 2.72 / 7 | 0.98 / 2 | 23 / 1 / 0 | 7.87 | 0.421 (0.418) | 0.355 / 0.353 | 0.529 | **79 209** |
| N2 | NEW | 2x8 | 0 | 2.72 / 7 | 0.98 / 2 | 23 / 1 / 0 | 7.87 | **0.565** (0.539) | 0.501 / 0.500 | 0.743 | 78 089 |
| B4 | BASE | 4x4 | 0 | 3.02 / 7 | 21.7 / 48 | 23 / 1 / 0 | 7.76 | 0.665 (0.645) | 0.591 / 0.565 | 0.773 | 48 961 |
| B2 | BASE | 2x8 | 0 | 3.02 / 7 | 21.7 / 48 | 23 / 1 / 0 | 7.76 | 1.013 (1.010) | 0.915 / 0.850 | 1.182 | 77 487 |

- `vgdrot-1011 vet_gd_halo_ipc: ON (node-local CUDA IPC, GPU-driven per-shell halo of the exact lists, no host wait)` is printed in N4, N4n, N2, C4g.
- Rotation cycles are no longer slower than steady ones in any healthy run (c10/c20 ~ 0.85-0.9x steady).
- Gas Newton fallbacks per cell-pass: C4 1.6e-2, C4g 2.7e-2; N4 4.7e-3; B4 4.3e-3. vimp positivity fallbacks: C4 41, C4g 31, others 0.
- Note: N4 and B4 both have time2_pred_cap limited cell-builds 1.094e7 of 3.86e8.

### conservative2 (C4 / C4g) vs matched (R4) and old conservative (RC4), bsgrec smoke 3359915
| | R4 matched | RC4 conservative | C4 conservative2 | C4g conservative2 + mg_gf |
|---|---|---|---|---|
| NON-CONVERGED | 0 | 0 | **6** | **4** |
| NOT ADMISSIBLE cycles | none | 2, 4-9 | 2-23 (all) | 2-23 (all) |
| first dt collapse | none | c3: 87.4 -> 0.0084 s | c3: 87.4 -> 0.0084 s, again c12 (0.0015), c19, c23 | c3: -> 0.0084 s, recovers to 15.8 s by c22, again c23 (0.0016) |
| min E from the solve | 7.76 | -9.19e3 | -4.05e6 | -2.97e4 |
| L_top after cycle 3 | 4.3-4.7e37 | **8.7e31** (lost) | 1.8-2.7e37 | see analysis.txt |
| t reached | 874 s (10 c) | 263 s (10 c) | 267 s (24 c) | see analysis.txt |
All four are identical through cycle 1 (L_top 1.67165e38). dt per cycle (C4): 87.37, 87.37, 87.37, 0.00836, 0.0167, 0.0334, 0.0669, 0.134,
0.268, 0.535, 1.07, 2.14, 0.00151, 0.00302, 0.00605, 0.0121, 0.0242, 0.0484, 0.0747, 0.0193, 0.0386, 0.0773, 0.116, 0.0184, 0.0368 s.
C4g: 87.37 x3, 0.00836, 0.0167, ..., 2.14 (c11), 0.515, 0.746, 1.08, 1.54, 2.21, 3.15, 4.47, 6.31, 8.89, 12.3, 15.8 (c22), 0.00162, 0.00323 s.
hst rows (t, dt, L_top, L_in, Etot, Picard) for C4 and C4g: deltaai-bsgaccel-1012/analysis.txt.

### GPU-h
4.2 min x 4 GH200 = **0.28 GPU-h**.

### Files
In git: docs/handover/deltaai-bsgaccel-1012/ (job.sub, analyse.sh, analysis.txt). On DeltaAI: /work/nvme/bivj/jma20/bsgaccel_1012/runs/3363179/<run>/.
