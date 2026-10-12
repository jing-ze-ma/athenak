# TASK raven -> DeltaAI: conservative2 + implicit_cs2_hll = drag on the full-size BSG and AG Car A (1 job, ≤ 1 h)

This task comes from Raven worker G (consfix-1011). It is a production gate: the user wants to know whether
`implicit_gas_src = conservative2` can be used on the BSG and on AG Car A.

- Why: on the reduced BSG (Raven CPU), conservative2 at 1fbe8621 had hesdirk2 stage 1 NOT ADMISSIBLE from
  cycle 2 on.
  - The bad cells were floor gas (rho 1e-20) with E/(rho c^2) ≈ 2-6.
  - There, the AP-HLL part of a face and the physical faces charged the gas all of d(aE)/dt, which gave kicks of
    0.2-0.5 c.
- The fix: commit **8cba89a6** (new key `implicit_cs2_hll = drag`, the default) charges only the drag on those
  faces.
  - Reduced BSG, 8 cycles: no fallback, dt the same as matched.
  - Reduced AG Car A: within 1e-7 of matched.

Files: docs/handover/consfix-c2drag-1012/ on bsg-files-1009.

## Limits
- 1 node x 4 GH200, 1 rank per GPU, `export MPICH_GPU_SUPPORT_ENABLED=1`, ≤ 1 h wall (job.sub stops new runs
  at 57 min).
- Reply SKIP if the queue wait is over about 2 h.

## 1. Build
- `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1).
- Source: fork branch `consfix-1011` @ **8cba89a6** (the coordinator pushes it).
- Report the md5.

## 2. Inputs (files in this directory; `@FILES@` → the resolved files directory)
- `bsgR.athinput.in`:
  - from bsg-files-1009 docs/handover/bsgrec-1011/bsgR.athinput.in;
  - @FILES@ = the BSG files used for the bsgrec smoke;
  - 7 lines added in `<rad_m1>`, all at defaults: implicit_cs2_hll, implicit_cs2_shift, dbg_t2_admiss,
    implicit_mg_levels, implicit_mg_halo, implicit_gf_modes, implicit_precond_float.
- `agcA.athinput.in`:
  - a copy of the current production deck /raven/ptmp/jinma/prodfinal_1011/run/agcA/agcar_A_fresh.athinput;
  - @FILES@ = docs/handover/agcarA_files (ic_agcar_A_ge.txt, rosseland_ext2_gs98_x0.36_z0.02.txt,
    planck_ext2orion_gs98_x0.36_z0.02.txt);
  - `implicit_cs2_hll` and `dbg_t2_admiss` added in `<rad_m1>`.
- Resolve both into `$G/bsgR.athinput` and `$G/agcA.athinput`, with G = /work/nvme/bivj/jma20/consfix_1011.
- `gpuwrap.sh`: copy the one from your bsgrec_1011 directory.

## 3. Runs (job.sub, in order; `DRY=1 BP=x bash job.sub` prints the commands)

| run | deck | cycles | keys |
|---|---|---|---|
| bM | BSG, 4 ranks x 4 blocks | 30 | implicit_gas_src = matched |
| bC2 | BSG | 30 | conservative2 (drag by default), dbg_t2_admiss = 2 |
| bC2g | BSG | 30 | bC2 + implicit_precond = mg_gf |
| aM | AG Car A, 4 ranks x 1 block | 60 | matched |
| aC2 | AG Car A | 60 | conservative2, dbg_t2_admiss = 2 |

All runs also get: `time/nlim`, `time/ndiag=1`, `output1/dt=1e-3` (hst every cycle) and no bin/rst dumps.

## 4. Report: NOTE-2026-10-12-deltaai-consfix-c2drag.md

**Per run:**
- rc, FATAL, NON-CONVERGED, the Picard line.
- The number of `hesdirk2 stage ... NOT ADMISSIBLE` lines (stage fallbacks) and their cycles.
- For any not-admissible stage, the `DBG_T2_ADMISS` blocks (min E / eint, r, rho, v).
- dt per cycle (from the `elapsed=` lines).
- s/cycle: steady cycles, excluding cycle 0, the first 2 cycles and multiples of 10.
- L_top per cycle (user.hst).

**BSG:**
- bC2 and bC2g vs bM: dt path.
- Max relative difference of L_top, E_rad and Etot at cycles 10, 20 and 30.

**AG Car A:** aC2 vs aM, per user.hst column, the max relative difference over the rows. Expect ≤ 1e-6.

**Verdict line:** conservative2 OK / not OK on BSG and on AG Car A.
