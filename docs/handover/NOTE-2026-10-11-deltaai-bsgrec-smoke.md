# NOTE DeltaAI: BSG true repro R / RC smoke + cycle-0 IC check (TASK-2026-10-11-deltaai-bsgrec-smoke) -- DONE

## RESULT (posted 2026-10-11 ~19:15 CDT; the job ran 09:05-09:07 CDT, the write-up was delayed on my side)
- **R (matched): clean.** R4 rc 0, FATAL 0, NaN 0, NON-CONV 0, Picard 3.95 / 7, no NOT ADMISSIBLE, steady 0.79 s/cycle (4 blocks/GPU).
- **RC (conservative): NOT healthy, do not start RC production as is.** rc 0 and NON-CONV 0, but hesdirk2 stage 1 is
  NOT ADMISSIBLE on every cycle 2-9 (min(E,T) < 0, the step is redone with backward Euler; 7 stage fallbacks, 8 BE steps of 10),
  the solve's min E is -9.19e3, vimp positivity fallbacks 6, Picard 15.4 / 51, and **dt collapses at cycle 3 from 87.37 s to
  0.0084 s** (then doubles back: 0.017, 0.019, 0.039, ... 0.62 s at cycle 10). **L_top falls from 4.46e37 to 8.7e31** (x2e-6).
  R and RC are identical to cycle 1 (same L_top 1.67165832e38), split at cycle 2.
- **Cycle-0 IC: BITWISE vs ic0_ref_1009.json, all 9 variables, for both R4 and RC4** (the original 10-09 IC).
- **Memory (hesdirk2, lowmem deck):** 4 blocks/GPU peak 48 961 MiB (R4) / 49 435 MiB (RC4) per GH200;
  **8 blocks/GPU (R2, Caltech layout) 77 497 MiB** -> fits a 141 GB H200 with room; 10-09 be/hr was 65.0 GiB.
- GPU-h: 1.3 min x 4 GH200 = **0.09 GPU-h**.


**Picked up** 2026-10-11 08:52 CDT.
- Built prodfinal-1011 **530c38e2** (= tip at build time): `athena_hes_gpu_530c38e2c428`, md5 **195a0a2fb73f509561b1e41df6b07014**.
- Files: our 10-09 files, md5 OK vs MD5SUMS; bsgR / bsgRC from the .in templates (diff = only implicit_gas_src + comments).
- ic0check.py sanity check on our own 10-09 GH200 hr cycle-0 dumps (bsg_hrdet_1009 job 3347471 detH_a, binary 98835d99):
  **BITWISE vs ic0_ref_1009.json in all 9 variables** (dens, velx, vely, velz, eint, m1_e, m1_f1, m1_f2, m1_f3).
- Job **3359915** (ghx4-interactive, 1 node x 4 GH200, 0:45): R4, RC4 (4 x 4, nlim 10), R2 (2 x 8, nlim 4), stop on rc != 0 / FATAL,
  nvidia-smi every 2 s. Extra command-line keys (both in the deck): `output1/dt=1e-3` (hst every cycle) and `output5/dt=1e30` (no rst).
  **PENDING** (priority 1987) at 09:00 CDT. STARTED line follows here when it runs.

(The STARTED line above was not pushed while the job ran; it ran 09:05:31-09:06:48 CDT on gh066.)

## Per run (job 3359915, gh066)
| run | layout | rc | FATAL / NaN / NON-CONV | Picard mean / max (solves) | dt | s/cycle steady (cycles 3-9, mean / median) | rotation cycle 10 | setup | GPU peak per GPU |
|---|---|---|---|---|---|---|---|---|---|
| R4 | 4 x 4 | 0 | 0 / 0 / 0 | 3.95 / 7 (19) | 87.37 s | 0.80 / 0.79 | 6.7 s (also the last cycle) | 5.2 s | 48 961 MiB |
| RC4 | 4 x 4 | 0 | 0 / 0 / 0 | 15.4 / 51 (19) | 87.37 -> **0.0084 s** at cycle 3 | 1.84 / 1.18 | 3.2 s | 4.3 s | 49 435 MiB |
| R2 | 2 x 8 (GPUs 0,1) | 0 | 0 / 0 / 0 | 5.0 / 7 (7) | 87.37 s | 1.48 (cycle 3 only) | - | 2.5 s | 77 497 / 77 487 MiB |

Solver lines:
- R4: `time_scheme=hesdirk2: stage steps=9 backward-Euler steps=1 stage fallbacks=0 vet clips=0 time2_vet_col=predict fallbacks=0
  time2_pred_cap=5.000000e-01 limited cell-builds=1.094451e+07 of 1.509949e+08`; vimp positivity fallbacks 0, min E 7.76.
- RC4: `stage steps=2 backward-Euler steps=8 stage fallbacks=7 ... time2_pred_cap=5.000000e-01 limited cell-builds=2.457600e+07 of 1.509949e+08`;
  `hesdirk2 stage 1 NOT ADMISSIBLE at cycle N (converged=1 min(E,T)=-1.000000e+00): the step is redone with backward Euler` for N = 2, 4, 5, 6, 7, 8, 9;
  vimp positivity fallbacks 6, min E from the solve -9.192631e+03; gas Newton fallbacks 5.5e7 (1.1e-2 per cell-pass; R4 1.2e-3).
- R2: stage steps 3, BE 1, fallbacks 0, limited cell-builds 1.094451e+07 of 5.033165e+07.
- No `implicit_gas_src` banner is printed by 530c38e2.
- The R4 rotation cycle (10) is also the last cycle, which holds the final output dumps, so 6.7 s is not a rotation cost.

## user.hst (every cycle; output1/dt=1e-3 on the command line)
R4: t, L_top, L_in, Etot, Picard
```
0.00000000e+00  0.00000000e+00  4.65121750e+37  1.56895735e+46  0.0
8.73700279e+01  1.67165832e+38  4.65121750e+37  1.56895630e+46  5.0
1.74739303e+02  4.45655061e+37  4.65121750e+37  1.56895584e+46  5.5
2.62108873e+02  4.32384983e+37  4.65121750e+37  1.56895588e+46  4.5
3.49478050e+02  4.74501134e+37  4.65121750e+37  1.56895588e+46  5.0
4.36846697e+02  4.55244711e+37  4.65121750e+37  1.56895589e+46  4.5
5.24214676e+02  4.68489328e+37  4.65121750e+37  1.56895589e+46  3.5
6.11581858e+02  4.63792351e+37  4.65121750e+37  1.56895589e+46  3.5
6.98948116e+02  4.64906632e+37  4.65121750e+37  1.56895589e+46  3.0
7.86313334e+02  4.64595730e+37  4.65121750e+37  1.56895589e+46  3.0
8.73677401e+02  4.64744008e+37  4.65121750e+37  1.56895589e+46  2.5
```
RC4: t, L_top, L_in, Etot, Picard
```
0.00000000e+00  0.00000000e+00  4.65121750e+37  1.56895735e+46  0.0
8.73700279e+01  1.67165832e+38  4.65121750e+37  1.56895630e+46  5.0
1.74739303e+02  4.45553525e+37  4.65121750e+37  1.56895584e+46  11.0
2.62108873e+02  3.04614651e+37  4.65121750e+37  1.56895638e+46  12.5
2.62117242e+02  6.99570316e+37  4.65121750e+37  1.56895638e+46  38.5
2.62133981e+02  1.10331657e+34  4.65121750e+37  1.56895634e+46  26.5
2.62153331e+02  6.51706207e+31  4.65121750e+37  1.56895634e+46  17.5
2.62192031e+02  8.69676527e+31  4.65121750e+37  1.56895631e+46  9.0
2.62269432e+02  8.69417400e+31  4.65121750e+37  1.56895631e+46  6.0
2.62424233e+02  8.68178989e+31  4.65121750e+37  1.56895631e+46  7.5
2.62733836e+02  8.66765106e+31  4.65121750e+37  1.56895630e+46  15.0
```
(L_top = 0 at t = 0 is the hst row written before the first step.)

## Files
In git: docs/handover/deltaai-bsgrec-1011/ (job.sub, analyse.sh, analysis.txt). On DeltaAI: /work/nvme/bivj/jma20/bsgrec_1011/runs/3359915/{R4,RC4,R2}
(out.txt, gpumem.log, user.hst, bin/ cycle-0 + final dumps).
