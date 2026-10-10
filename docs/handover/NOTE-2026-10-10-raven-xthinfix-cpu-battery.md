# NOTE Raven CPU battery: xthinfix-1009 e50ec7d6 new arms (sqs, sqb, sqp, sqsb, sqsbp, aqs), 10-10

This NOTE answers TASKS.md (the Raven takeover from viper). It also covers the plm-gate rerun from TASK-2026-10-10-caltech-plm-gates-rerun.
- Source: /raven/ptmp/jinma/xthinfix_1010/src (fork/xthinfix-1009 e50ec7d6, kokkos 4.6.2).
- Run tree: /raven/ptmp/jinma/xthinfix_1010/cpu/run (runs/, gates/, beams/, RESULTS/).
- Scripts: /raven/ptmp/jinma/xthinfix_1010/cpu (battery_raven.sh, run.sbatch, fixkeys.py, ana/ with the od.py copy).

STATUS: complete. Jobs (Raven, partition general, 72 cores per node, OMP 1 per run):
- Build: 31050045 (interactive partition, 1.6 min).
- Smoke: 31050147 (beam ba0 sqsb) and 31050148 (gates sqsb), interactive partition. No key FATAL.
- Order: 8 round-robin parts, 31050182-31050193. All 1550 runs finished with rc 0 in about 11 min.
- Gates: 31050195 (88 runs). Beams: 31050197 (50 runs). Ana: 31050563.
- The held agcA_rc / agcB_rc jobs were not touched.

## Binary

- `bin/athena_cpu_e50ec7d6`: built-in pgens, MPI (gcc/13, openmpi/5.0, cmake/3.30), Release, Kokkos OpenMP host, no CUDA.
- md5 `48da95b9cee7904eeadfc6f60943aefa`.
- Reproduces Caltech 413b38af: every reference-arm number below that Caltech also ran (cen / hr / sq / aq order tables, sq / aq beams) is digit-identical to NOTE-2026-10-10-caltech-arm-matrix. That includes sq ba0 diverging at cycle 4 and ba20 at cycle 88, and xb20 sq Picard 65.2 (max 177).

## Verdict

0. **Acceptance criteria** (user, TASKS.md):
   - (1) Second order in X / T / C. X order of about 1 is accepted only in the thin part (tau < 0.1 or Kn > Kn0); orders of about 0, negative orders, or errors that grow with resolution FAIL.
   - (2) Passes the simple tests (pulses, rw, atm vs Hopf, G1 / G3 / G5) and the beams, with no divergence and no NON-CONVERGED.
   - (3) Cheap: Picard close to hr.

   **No new arm passes. Each new arm fails (2) on the beams and (3) on cost, and the sq family also inherits sq's failure of (1) on rw_t10.**

1. **The new keys do almost nothing in 1-D order and in the gates.**
   - Pulses (all kappa) and rw_t1000: sqs, sqb, sqp, sqsb and sqsbp are bitwise equal to sq (orders and L1).
   - rw_t10: the new arms differ only in the third digit.
   - atm: sqb / sqp / sqsbp equal sq. sqs / sqsb change the thin-region self-convergence (item 3).
   - aqs equals aq except in the atm thin top.
   - So the sq family keeps sq's convergence:
     - thin pulse k0.128 hesdirk2 X 1.75 1.94 2.04 (= cen);
     - k12.8 X 2.25 2.24; thick k1280 X 1.93;
     - atm thick (tau > 10) 1.60, tau 1-10 1.15, tau < 0.1 1.12, Kn > 0.3 0.95;
     - Hopf L1 at 512 is 3.87e-4 (sq / sqb / sqp / sqsbp) and 3.79e-4 (sqs / sqsb).
   - aqs keeps aq's failures:
     - thin pulse k0.128 be X 1.76 1.76 **-0.83**, hesdirk2 X 1.74 1.91 **0.68**;
     - k12.8 be X **-0.58**;
     - rw_t10 E hesdirk2 X 1.11 **-1.47 -1.00**, C **-0.96 -0.47**, be X **-0.68 -1.61 -0.83**.

2. **The sq family fails (1) on rw_t10 E, like sq itself.** This was not flagged at Caltech, which quoted only the hesdirk2 X orders.
   - be X: 2.17 2.0x **-3.19**. L1 at 512 is 1.14e-2, against cen 3.17e-4. The error grows 9x from 256 to 512, for all six sq arms.
   - hesdirk2 T: 1.00 **-0.12 to -0.16** (cen 1.01 1.00 1.00). The error stalls at 5.2e-4.
   - be C: last order 0.58 (cen 0.87).
   - hesdirk2 X (1.76-1.97) and dens converge.
   - Among the new arms, sqs / sqsb / sqb have slightly better hesdirk2 X last orders (1.94-1.97 vs sq 1.88). sqp / sqsbp are slightly worse (1.76-1.79).

3. **step mode (sqs, sqsb) makes the atm thin regions non-monotone.**
   - tau 0.1-1: e 5.42e-3, 3.51e-4, **7.10e-4**, 4.23e-4. Orders 3.95 **-1.01** 0.75; the error grows from 64 to 128.
   - Kn > 0.3: orders 2.41 **0.24** 0.97.
   - The finest errors are comparable to sq's (4.2e-4 vs 4.9e-4).
   - Under criterion (1) this is a marginal FAIL, because the error grows at one level. aqs shows a smaller effect (tau < 0.1: 0.37 1.61 1.17).

4. **Gates: everything passes, and the plm segfault is fixed.**
   - All 88 gate runs have rc 0 and FATAL 0. The plm-gate rerun (sq, aq, knp on G1 / G3, which were rc 139 at Caltech) is now rc 0.
   - Viper's sq check is reproduced: G1 cfl 1e2 max|dE/E| 1.738e-2, L1 6.514e-4 (viper 1.74e-2 / 6.5e-4); G3 cfl 10 L1 0.0138 (viper 0.0138).
   - G1 (all arms except cen): max|dE/E| 1.63e-2 to 1.77e-2.
     - L1 at 1e2 / 1e4: sq family 6.5e-4 / 7.2e-4; aq / aqs 1.21e-3 / 1.94e-3; knp 7.8e-4 / 1.24e-3; hr 1.56e-3 / 1.95e-3.
     - cen gives max|dE/E| 0.17 (known: the central flux in the thin top).
   - G3 L1(E) 0.0135-0.0140: PASS for every arm.
   - G5 k1280 / k12800: PASS for every arm (ratio 1.0167 / 0.99996).
     - k12.8 / k128 FAIL (ratio 0.53-0.71) for every arm including cen and hr, as at Caltech. These are not in the diffusion limit.
   - NON-CONVERGED in G1 cfl 1e4 (gate input implicit_maxit = 30): every arm has 25-33 non-converged steps (cen 25, hr / aq / knp 32, sq family 33).
     - sqp / sqsbp: 44-45, with Picard mean 29.5 out of a cap of 30. The predictor adds one dc phase per step: +1.0 Picard in G1 1e2 / G3 / pulses, and +7.5 at G1 1e4.
     - G1 1e2, G3 and G5 have NON-CONVERGED 0, apart from 1 for cen in G1 1e2.

5. **Beams (100 cycles, c dt/dx ~ 300): the new sq arms are LESS stable than sq. ba0 and ba20 still diverge.** Each divergence ends as a FATAL "Picard solve DIVERGED".
   - sq (reference): ba0 diverges at cycle 4, ba20 at cycle 88; xb20 Picard 65.2 (max 177).
   - sqs: diverges on xb20 (c8), ba0 (c5) and ba20 (c11). Survives cyl and shd3b at 13.5 / 10.5 Picard (sq 9.8 / 4.7).
   - sqb: the most stable sq variant.
     - Survives ba20 (Picard 39.5) and xb20 (43.6, vs sq 65.2).
     - Still diverges on ba0, but late (cycle 74 instead of 4).
     - cyl 18.6, shd3b 6.3.
   - sqp: diverges on xb20 (c44), ba0 (c5) and ba20 (c12). cyl 21.9, shd3b 14.4.
   - sqsb: diverges on 4 of 5 (xb20 / ba0 / ba20 at c3, shd3b at c4). Only cyl survives (30.4).
   - sqsbp: diverges on 4 of 5 (c2-c4). Only cyl survives (34.0, max 183).
   - aqs: stable on all 5 (as aq), but more expensive than aq: xb20 12.3 / ba0 17.7 / ba20 26.5 / cyl 33.8 / shd3b 7.5, against aq 7.8 / 15.4 / 22.4 / 33.8 / 4.5.
   - hr: 4.0 / 4.7 / 5.0 / 5.0 / 3.0 with BiCGStab inner mean 17-45.
   - Every run that finished has NON-CONVERGED 0.

6. **Cost: the step mode (qs_mode = step) does NOT lower Picard; it raises it.**
   - Picard, step vs pass: shd3b sqs 10.5 vs sq 4.7; cyl 13.5 vs 9.8; aqs vs aq +1 to +4 per beam; atm 1-D 5.7 vs 4.2.
   - pred adds about 1 Picard pass per step everywhere, and up to +10 on the beams.
   - BiCGStab inner mean per Picard pass:
     - plm arms 5-21, about the same as sq;
     - hr 17-45; cen 47-106 (max 200).
   - Per step the plm arms remain 2-10x hr in Picard (cyl 13-34 vs 5) wherever they survive.
   - The 1-D order runs are not discriminating: Picard mean of means is 1.9-2.1 for all arms, 3.1 with pred.
   - Wall times come from shared-node runs (50-72 concurrent jobs), so they are not timings.

7. **FATAL / NON-CONVERGED summary**
   - Order (1550 runs): FATAL 0, NON-CONVERGED 0.
   - Gates (88 runs): FATAL 0. NON-CONVERGED only in G1 cfl 1e4 (every arm, maxit 30) and 1 in cen G1 1e2.
   - Beams: 17 FATAL "Picard solve DIVERGED", all in sq-family arms: sq 2, sqs 3, sqb 1, sqp 3, sqsb 4, sqsbp 4. The smoke tree has 1 more, a duplicate of ba0 sqsb.
   - No key FATAL anywhere.

### PASS / FAIL against the acceptance criteria

Per-region X orders are the last two levels.
- The thin region is pulse_k0.128 hesdirk2 X plus atm tau < 0.1 / Kn > 0.3.
- The intermediate region is pulse_k12.8 hesdirk2 X plus atm tau 0.1-1 / 1-10.
- The thick region is pulse_k1280 X plus atm tau > 10.

| arm | thin X | intermediate X | thick X | rw_t10 E (be X / hd2 X / hd2 T last 2) | T vs cen | gates G1/G3/G5thick | beams | Picard vs hr | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cen | 1.94 2.04 / atm 1.01, 0.94 | 2.25 2.24 / 0.97, 0.81 | 1.93 / 1.91 | 2.04 1.96 / 1.97 1.90 / 1.00 1.00 | ref | G1 FAIL (Hopf 0.17), G3 / G5 PASS | 5/5 stable | 0.5-1x, inner 2-4x | reference (fails thin Hopf) |
| hr | **0.69** / atm 1.70, **-0.01** | 1.07 / **0.27, 0.06** | 1.81 / **-0.42** | **-1.53 -0.79** / **-1.67 -1.18** / **-0.07 0.34** | worse | PASS | 5/5 stable | 1x | FAIL (1) |
| sq | 1.94 2.04 / 1.12, 0.95 | 2.24 / 0.84, 1.15 | 1.93 / 1.60 | 2.02 **-3.19** / 1.97 1.88 / **-0.15 -0.13** | = cen except rw_t10 | PASS | **ba0 c4, ba20 c88 DIV** | 1.5-16x | FAIL (1) rw_t10, (2), (3) |
| aq | **0.68** / 1.09, 1.00 | 1.07 / 0.97, 1.11 | 1.93 / 1.58 | **-1.61 -0.83** / **-1.47 -1.00** / 0.02 0.37 | worse | PASS | 5/5 stable | 1.5-7x | FAIL (1) |
| sqs | = sq / 1.14, 0.97 (mid 0.24) | = sq / **-1.01** 0.75, 1.13 | = sq / 1.60 | 2.01 **-3.18** / 1.97 1.97 / **-0.12 -0.16** | as sq | PASS | **xb20 c8, ba0 c5, ba20 c11 DIV** | 2.2-3.5x | FAIL (1), (2), (3) |
| sqb | = sq | = sq | = sq | 2.01 **-3.19** / 1.97 1.94 / **-0.16 -0.13** | as sq | PASS | **ba0 c74 DIV** (ba20, xb20 survive) | 2-11x | FAIL (1), (2), (3) |
| sqp | = sq | = sq | = sq | 2.03 **-3.20** / 1.97 1.79 / **-0.12 -0.14** | as sq | PASS (G1 1e4 NC 45) | **xb20 c44, ba0 c5, ba20 c12 DIV** | 4.4-5x | FAIL (1), (2), (3) |
| sqsb | as sqs | as sqs | = sq | 2.01 **-3.18** / 1.97 1.96 / **-0.12 -0.16** | as sq | PASS | **4 of 5 DIV (c3-c4)** | 6x (cyl) | FAIL (1), (2), (3) |
| sqsbp | = sq | = sq | = sq | 2.03 **-3.20** / 1.97 1.76 / **-0.12 -0.16** | as sq | PASS (G1 1e4 NC 44) | **4 of 5 DIV (c2-c4)** | 7x (cyl) | FAIL (1), (2), (3) |
| aqs | **0.68** / 1.17, 1.05 | 1.07 / 0.99, 1.14 | 1.93 / 1.59 | **-1.61 -0.83** / **-1.47 -1.00** / 0.02 0.37 | worse | PASS | 5/5 stable | 2.5-7x | FAIL (1), (3) |

Notes on the table:
- The atm entries are the last two orders (128 to 256, 256 to 512) of ana_atm_self.py.
- Pulse entries are hesdirk2 X, except thick, which is be X / atm.
- **T order of about 2 is not reached by cen either in this battery.**
  - hesdirk2 T is 1.0 for cen on rw_t10, rw_t1000 and pulse_k12.8.
  - It is 1.57-1.79 on pulse_k0.128 and 1.90-1.96 on pulse_k1280.
  - be is first order by construction.
  - So the T criterion is judged relative to cen: the sq family equals cen except on rw_t10 E.

## Order tables

From RESULTS/order_eval.txt: E orders, with dens for rw. `= sq` means bitwise equal to sq. Full raw output: /raven/ptmp/jinma/xthinfix_1010/cpu/run/RESULTS/order_eval.txt. The compact form comes from cpu/ana/ordtab.py.

| case sch mode | cen | hr | sq | aq | sqs | sqb | sqp | sqsb | sqsbp | aqs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| k0.128 be X | 1.76 1.97 2.15 | 1.76 1.76 -0.82 | 1.76 1.97 2.15 | 1.76 1.76 -0.83 | = sq | = sq | = sq | = sq | = sq | = aq |
| k0.128 hd2 X | 1.75 1.94 2.04 | 1.74 1.90 0.69 | 1.75 1.94 2.04 | 1.74 1.91 0.68 | = sq | = sq | = sq | = sq | = sq | = aq |
| k0.128 hd2 T | 1.79 1.71 1.57 | 1.24 1.13 1.06 | = cen | 1.23 1.12 1.06 | = sq | = sq | = sq | = sq | = sq | = aq |
| k0.128 hd2 C | 1.75 1.95 2.04 | 1.75 1.93 1.72 | = cen | 1.75 1.93 1.72 | = sq | = sq | = sq | = sq | = sq | = aq |
| k12.8 be X | 1.89 2.25 2.25 | 1.89 2.18 -0.55 | = cen | 1.89 2.19 -0.58 | = sq | = sq | = sq | = sq | = sq | = aq |
| k12.8 hd2 X | 1.90 2.25 2.24 | 1.90 2.24 1.07 | = cen | 1.90 2.24 1.07 | = sq | = sq | = sq | = sq | = sq | = aq |
| k12.8 hd2 T | 1.05 1.03 1.02 | 1.03 1.02 1.01 | = cen | 1.03 1.02 1.01 | = sq | = sq | = sq | = sq | = sq | 1.03 1.01 1.01 |
| k1280 be/hd2 X | 1.89 1.97 1.93 | 1.89 1.95 1.80 | = cen | = cen | = sq | = sq | = sq | = sq | = sq | = cen |
| k1280 hd2 T | 1.96 1.94 1.90 | 2.08 2.55 -0.97 | = cen | 1.96 1.95 1.93 | = sq | = sq | = sq | = sq | = sq | 1.98 1.98 1.97 |
| rw_t10 E be X | 2.17 2.04 1.96 | -0.79 -1.53 -0.79 | 2.17 2.02 -3.19 | -0.68 -1.61 -0.83 | 2.17 2.01 -3.18 | 2.17 2.01 -3.19 | 2.17 2.03 -3.20 | 2.17 2.01 -3.18 | 2.17 2.03 -3.20 | = aq |
| rw_t10 dens be X | 2.18 2.14 2.16 | 2.09 1.22 -0.03 | 2.18 2.14 1.62 | 2.10 1.23 -0.06 | = sq | = sq | 2.18 2.13 1.63 | = sq | 2.18 2.13 1.63 | = aq |
| rw_t10 E hd2 X | 1.91 1.97 1.90 | 1.54 -1.67 -1.18 | 1.91 1.97 1.88 | 1.11 -1.47 -1.00 | 1.91 1.97 1.97 | 1.91 1.97 1.94 | 1.91 1.97 1.79 | 1.91 1.97 1.96 | 1.91 1.97 1.76 | = aq |
| rw_t10 E hd2 T | 1.01 1.00 1.00 | -0.86 -0.07 0.34 | 1.00 -0.15 -0.13 | -0.49 0.02 0.37 | 1.01 -0.12 -0.16 | 1.01 -0.16 -0.13 | 1.00 -0.12 -0.14 | 1.00 -0.12 -0.16 | 0.99 -0.12 -0.16 | = aq |
| rw_t10 E hd2 C | 2.35 2.42 1.12 | 0.73 -1.22 -0.65 | 2.35 2.41 0.74 | 0.42 -0.96 -0.47 | 2.35 2.41 0.72 | 2.35 2.43 0.72 | 2.35 2.42 0.71 | 2.35 2.40 0.75 | 2.35 2.40 0.74 | = aq |
| rw_t10 E be C | 0.59 0.74 0.87 | 0.00 0.31 2.67 | 0.56 0.79 0.58 | 0.01 0.31 2.69 | = sq | = sq | 0.59 0.74 0.58 | 0.57 0.78 0.58 | 0.59 0.74 0.58 | = aq |
| rw_t1000 E X (both) | 2.06 1.92 | 1.96 1.34 | = cen | = cen | = cen | = cen | = cen | = cen | = cen | = cen |
| rw_t1000 T / C | 0.97-1.07 | same | = cen | = cen | = cen | = cen | = cen | = cen | = cen | = cen |

atm (RESULTS/atm_hopf.txt, atm_self_regions.txt). Self-convergence orders run 32 to 512. Hopf L1 and top-cell |dE/E| are at 512.

| arm | Hopf L1 | top | self E | tau<0.1 | 0.1-1 | 1-10 | >10 | Kn>0.3 | Picard (1-D atm) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cen | 8.40e-3 | 1.5e-1 | 1.99 1.96 1.61 | 1.06 1.02 1.01 | 0.95 0.87 0.97 | 2.39 1.79 0.81 | 2.06 2.13 1.91 | 0.86 0.92 0.94 | 2.7 |
| hr | 4.24e-3 | 8.4e-4 | 1.10 0.13 -0.29 | 0.26 0.79 1.70 | 0.12 -0.34 0.27 | 0.42 0.06 0.06 | 1.26 0.17 -0.42 | 0.14 0.12 -0.01 | 2.4 |
| sq / sqb / sqp / sqsbp | 3.87e-4 | 4.1e-6 | 1.90 1.76 1.41 | 1.08 1.11 1.12 | 1.34 0.80 0.84 | 1.75 1.44 1.15 | 1.98 1.94 1.60 | 1.34 0.97 0.95 | 4.2 / 3.8 / 6.9 / 6.8 |
| sqs / sqsb | 3.79e-4 | 4.6e-5 | 2.12 1.60 1.40 | 1.52 0.70 1.14 | 3.95 **-1.01** 0.75 | 2.40 1.11 1.13 | 2.07 1.90 1.60 | 2.41 **0.24** 0.97 | 5.7 / 4.4 |
| aq | 3.89e-4 | 4.8e-6 | 1.87 1.74 1.40 | 1.14 1.12 1.09 | 1.26 0.97 0.97 | 1.67 1.40 1.11 | 1.96 1.93 1.58 | 1.29 1.05 1.00 | 4.3 |
| aqs | 3.81e-4 | 4.5e-5 | 1.83 1.79 1.42 | **0.37** 1.61 1.17 | 1.16 1.14 0.99 | 1.63 1.48 1.14 | 1.96 1.94 1.59 | 1.00 1.31 1.05 | 4.9 |

(The sqp / sqsbp atm numbers differ from sq only in the third digit.)

## Gates (RESULTS/gates.txt)

G1 cells give Picard mean, NON-CONVERGED, max|dE/E| and L1 vs Hopf. G3 cells give L1(E) and Picard. G5 lists the ratio for k12.8 / 128 / 1280 / 12800; the target is 1 +- 0.02.

| arm | G1 cfl 1e2 | G1 cfl 1e4 | G3 c1 | G3 c10 | G5 ratio |
| --- | --- | --- | --- | --- | --- |
| cen | 2.1, NC 1, 1.70e-1, 8.66e-3 | 16.7, NC 25, 1.70e-1, 8.66e-3 | PASS 0.0135, 2.01 | PASS 0.0136, 2.53 | 0.715 0.534 1.0167 0.99996 |
| hr | 2.1, NC 0, 1.75e-2, 1.56e-3 | 21.5, NC 32, 1.77e-2, 1.95e-3 | PASS 0.0138 | PASS 0.0140 | 0.566 0.529 1.0166 0.99997 |
| sq | 2.1, NC 0, 1.74e-2, **6.51e-4** | 22.0, NC 33, 1.63e-2, 7.24e-4 | PASS 0.0136 | PASS **0.0138** | 0.552 0.534 1.0167 0.99996 |
| aq | 2.1, NC 0, 1.75e-2, 1.21e-3 | 21.5, NC 32, 1.77e-2, 1.94e-3 | PASS 0.0137 | PASS 0.0140 | 0.553 0.529 1.0166 0.99996 |
| knp | 2.1, NC 0, 1.74e-2, 7.82e-4 | 21.5, NC 32, 1.77e-2, 1.24e-3 | PASS 0.0136 | PASS 0.0138 | 0.552 0.534 1.0167 0.99996 |
| sqs | = sq | = sq | PASS 0.0136 | PASS 0.0138 | 0.552 0.534 1.0167 0.99996 |
| sqb | = sq | = sq | PASS 0.0136 | PASS 0.0138 | 0.552 0.534 1.0167 0.99996 |
| sqp | 3.1, NC 0, = sq | 29.5, NC 45, 1.63e-2, 7.24e-4 | PASS 0.0136, 3.01 | PASS 0.0138, 3.53 | = sq |
| sqsb | = sq | = sq | PASS 0.0136 | PASS 0.0138 | = sq |
| sqsbp | 3.1, NC 0, = sq | 29.4, NC 44, = sqp | PASS 0.0136, 3.01 | PASS 0.0138, 3.53 | = sq |
| aqs | = aq | = aq | PASS 0.0137 | PASS 0.0140 | = aq |

All gate runs have rc 0. The G3 Picard means are 2.01 (c1) and 2.53 (c10) for every arm without pred.

## Beams (RESULTS/beams.txt)

100 cycles, np 1, c_light 1000. Each cell gives Picard mean / max, BiCGStab inner mean, and wall time (shared node, not a timing). DIV is a FATAL "Picard solve DIVERGED" at the given cycle.

| arm | xb20 | ba0 | ba20 | cyl | shd3b |
| --- | --- | --- | --- | --- | --- |
| cen | 3.7/5, 95, 578s | 4.7/5, 106, 757s | 4.2/5, 80, 564s | 2.2/6, 47, 138s | 3.0/4, 75, 633s |
| hr | 4.0/4, 37, 380s | 4.7/5, 40, 430s | 5.0/6, 45, 493s | 5.0/5, 17, 141s | 3.0/3, 21, 363s |
| sq | 65.2/177, 7, 1187s | **DIV cycle 4** | **DIV cycle 88** | 9.8/74, 7, 153s | 4.7/14, 18, 562s |
| aq | 7.8/12, 15, 368s | 15.4/24, 17, 624s | 22.4/29, 17, 852s | 33.8/40, 6, 318s | 4.5/14, 17, 426s |
| sqs | **DIV cycle 8** | **DIV cycle 5** | **DIV cycle 11** | 13.5/38, 7, 182s | 10.5/19, 17, 957s |
| sqb | 43.6/77, 6, 758s | **DIV cycle 74** | 39.5/70, 10, 943s | 18.6/89, 6, 213s | 6.3/16, 12, 551s |
| sqp | **DIV cycle 44** | **DIV cycle 5** | **DIV cycle 12** | 21.9/41, 9, 272s | 14.4/22, 21, 1431s |
| sqsb | **DIV cycle 3** | **DIV cycle 3** | **DIV cycle 3** | 30.4/82, 5, 284s | **DIV cycle 4** |
| sqsbp | **DIV cycle 2** | **DIV cycle 2** | **DIV cycle 3** | 34.0/183, 6, 332s | **DIV cycle 4** |
| aqs | 12.3/24, 12, 450s | 17.7/26, 17, 688s | 26.5/41, 17, 977s | 33.8/40, 7, 343s | 7.5/15, 19, 571s |

Other beam diagnostics:
- plm positivity fallbacks (finished runs):
  - xb20: sq 6512, sqb 4346, aq 736, aqs 1186.
  - ba20: aq 2230, aqs 2620, sqb 3918.
  - shd3b: sqs 973, sqp 1051 (sq 373).
- sqb prints "implicit_hr_pos = bound: cell-axes set to dc=..." and has more `solved E <= e_floor` cell-solves than sq on cyl (1321 vs 0).
- sqp prints "implicit_hr_recon_pred: dc-phase passes / plm-phase passes" (shd3b 296 / 1145).
- Accuracy (cylref / collref / shdfs) of the survivors is close to sq / aq.
  - cyl <|E/J-1|> 0.37-0.43. sqs reaches face max |F|/cE 3.21 and sqp 1.72; the others are at most 0.65.
  - xb20 beam sum E / sum J: 2.78-2.90 (hr 2.82).
  - shd3b E/ex lit: sqs 1.515 and sqp 0.36 at x 0.65 vs sq 0.556; this is noisier.
  - Full output is in beams.txt.

## Raven adaptations (copies in cpu/, xtf/ and src/ untouched)

- **fixkeys.py**: a command-line key missing from a fresh input is FATAL.
  - Every gate / beam run gets `<run>/in.athinput`, a copy of its input with every missing command-line key added. The list goes to fixkeys.txt.
  - Additions: vet_tensor in atm2d, and implicit_hr_recon_qs_mode / implicit_hr_recon_pred everywhere they were passed.
  - The generated pulse2d / marsh2d inputs also carry `vet_tensor = full`, `implicit_hr_recon_qs_mode = pass` and `implicit_hr_recon_pred = false` in the gate key block.
  - The order inputs need nothing: od.py write_input appends every rad key to `<run>/in.athinput`.
- **od.py copy** (cpu/ana/od.py): ARMS gains sqs, sqb, sqp, sqsb, sqsbp and aqs, with the same keys as battery_raven.sh.
  - Runner semantics, checked: `run` reverses the job file (finest first), skips any run dir that already has wall.txt (so a job can be resumed), and runs one.sh with OMP 1 under nice.
  - Order was split into 8 round-robin files (awk NR%8), one general node each with P = 72.
- **/usr/bin/time** is replaced by a `date +%s.%N` stand-in that writes the same "%e s wall" line to time.txt.
- **cmd.txt**: the gate / beam worker writes `<run>/cmd.txt` (binary, input and keys) for cylref / collref.
- **python**: the system python3 and the miniconda base have no numpy. `cpu/pybin/python3` links to /u/jinma/miniconda3/envs/magritte/bin/python (numpy 1.23.5, scipy, h5py) and comes first in PATH.
- **gates.out** lines now also carry fatal, wall, Picard max / NON-CONVERGED and inner iterations.
- **Arms**:
  - gates: cen hr sq aq knp + 6 new;
  - beams: cen hr sq aq + 6 new;
  - order: cen hr sq aq + 6 new.
  - kn01 / kn03 / kn1 / hrx0 / sqf / sqv / sqa were not rerun.
- Smoke: smoke/ (ba0 sqsb diverged at cycle 3, the same result as the full run; gates sqsb all rc 0).
