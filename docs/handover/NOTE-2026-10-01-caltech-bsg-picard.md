# NOTE 2026-10-01 (Caltech): BSG port gate -- same Picard-mean flip, here between CPU and GPU on one machine

Adds a Caltech data point to NOTE-2026-10-02-deltaai-bsg-picard-question / NOTE-2026-10-02-viper-bsg-picard-answer.

## Setup
- Built from rt-integration 8b2d770c (= 30bf6c03 + handover docs), he_star_m1. CPU: x86 gcc 13.2 + HPC-X OpenMPI.
  GPU: H200 CUDA 12.9 sm_90.
- Column bsg_col_arm2 from the bundle, 1 rank. c1 = CPU, g1/g2 = GPU (gate job 3748257).

## Result
| run | Picard mean | dt_end | L_top/L_in end | FATAL/nan |
|---|---|---|---|---|
| c1 (x86 CPU) | **2.994** (= viper) | 1.2430144e+02 | 0.999466 | 0/0 |
| g1, g2 (H200) | **3.000** (= DeltaAI CPU) | 1.2430144e+02 | 0.999466 | 0/0 |

- c1: one step takes 2 passes. It is in hst row 10 (t 9000-10000 s, row mean 2.875 = 1 of 8 steps), so it is
  consistent with viper's step 78 (t ~ 9.7e3 s). Our x86 CPU reproduces viper's x86 CPU.
- c1 vs g1, despite the different pass count: hst dt/L_in/E_rad/e_gas <= 7e-15, L_bot 3e-14, L_top 1e-10;
  dumps dens/eint <= 4e-19, velx (rho > 1e-17) 3.7e-8, m1_e 0. g1 vs g2 bitwise (hst + all dumps).
- So the flip is not specific to x86 vs Grace: it also happens between CPU and GPU on one machine. That fits viper's
  reading of step 78 (a lucky first linear solve, bound 1.5e-10).

## Gate change at Caltech
Our automatic checker had a CPU-vs-GPU Picard tolerance of 5e-4, so the gate FAILED on this alone. With the user's
OK it is now 0.02 (one flipped pass), in line with viper's "Picard mean report-only". Other criteria are unchanged
(dt_end, L_top/L_in, 0 FATAL/nan, CPU vs GPU <= 1e-5 hst / dumps). The re-check gives GATE_PASS. The 3-D smoke
(20 cycles, 4 ranks) runs next in gate job 3750607 (2 nodes x 2 H200); on pass it touches READY for production
3746142 (4 nodes x 1 H200). Status NOTE-2026-10-0x-caltech-bsg.md follows.

**Other sites:** if your gate compares the Picard mean between CPU and GPU with a tight tolerance, expect the same
flip. Treat it as report-only.
