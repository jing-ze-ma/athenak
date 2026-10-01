# NOTE 2026-10-01 (DeltaAI): BSG arm 2 port gate PASSED; production NOT run (user decision)

Reply to TASK-2026-10-01-deltaai-bsg.md. Machine: NCSA DeltaAI, 4x GH200 per node, CUDA.

**Status:** the port gate (4a + 4b) passed. **The user dropped the DeltaAI production on 10-01 ~12:15 CDT**: they want
the BSG production run in one go, not in pieces. On DeltaAI a 1-node 24 h ghx4 job had a `--test-only` start
estimate of **2026-10-11**. Only short `ghx4-interactive` jobs start the same day (2 h cap, 1 job per user, and
in practice only <= 45 min fit the backfill). No production link ran. Nothing is queued on DeltaAI.

## 1. Binaries (built before NOTE-2026-10-01-bsg-arm2-in-rt-integration, so from bsg-arm2, NOT rt-integration)

| | commit | md5 |
|---|---|---|
| GPU (CUDA, sm_90 + Grace, MPI) | bsg-arm2 d0abe21e | 44b69397a61bf49b77e524bc971c5f81 |
| CPU (host-only Grace, no MPI; port gate) | bsg-arm2 d0abe21e | e11f1f9fc56931fa4894e42322f14d3a |

- **Stack:** the bench-2026-09-29-hebox stack, via `build_inc_deltaai.sh` with new targets `hes_gpu` / `hes_cpu`.
- **No port fix needed:** the c2p_track inline-lambda trap did not bite, since the BSG run is hydro. vet_col on CUDA
  gave no scratch-size error.
- **Bundle issue:** `bsg_1001_bundle/SETUP.sh` runs `md5sum -c MD5SUMS` relative to the current directory, so it
  fails unless called from inside the bundle directory.

## 2. Gate 4a: 1-D column (`bsg_col_arm2.athinput`, 161 cycles to 2e4 s)

`gate_compare.py` summary lines (identical for c1 = CPU on the login node, and g1, g2 = 1 GH200 each):

```
rows 21  t_end 2.000000e+04  dt_end 1.243014e+02  Picard mean 3.000  L_top/L_in end 0.999466 range [0.999373, 0.999467]  FATAL 0  nan 0
```

These match viper's CPU and MI300A numbers to every printed digit. Newton fallbacks: 0 on CPU.

| c1 vs g1 (max rel diff) | DeltaAI | viper CPU vs MI300A |
|---|---|---|
| hst dt / L_top / L_bot | 1.3e-8 / 4.4e-10 / 2.1e-7 | 1.0e-8 / 6.8e-10 / 2.6e-7 |
| hst E_rad / e_gas | 1.4e-8 / 1.3e-8 | 1.7e-8 / 1.6e-8 |
| hst KE_int / v1sq_wall / Min_top | 8.7e-6 / 8.9e-4 / 2.9e-4 | 1.9e-5 / 1.8e-3 / 2.3e-4 |
| dump dens / eint / m1_e / m1_f1 | 1.7e-8 / 9.2e-8 / 1.2e-7 / 3.8e-5 | 6.8e-8 / 9.2e-8 / 1.2e-7 / 3.8e-5 |
| dump velx: all cells / rho > 1e-17 | 7.3e-2 / 3.9e-6 | 5.8e-2 / 4.0e-6 |

**g1 vs g2: bitwise identical** in every hst column. Wall time is about 30 s per GPU run, including startup. **PASS.**

## 3. Gate 4b: 3-D smoke (`bsg3d_arm2.athinput`, 1 node x 4 GH200, 1 block per GPU, 20 cycles)

- **Job and exit:** job 3286515, rc 0, 0 FATAL/nan.
- **Time step:** dt 88.1765 s (viper 88.18).
- **Solver:** Picard 3.0 mean and 3 max, NON-CONVERGED 0, Newton fallbacks 0.
- **Speed: median 0.67 s/cycle** over cycles 5-20 (about 0.55 s/cycle steady from cycle 9). Viper's MI300A
  timing job gave 1.15 s/cycle, so GH200 is about 1.7-2x faster per GPU at the production load.
- **Projection:** a full arm-2 run of about 56,300 cycles would take about 9-11 h on 1 DeltaAI node, if it could
  get a node.

**PASS.**

## 4. Production

Not run (user, see top). The guarded link and chain scripts are ready in `/work/nvme/bivj/jma20/bsg_1001`
(`prod.sub`, `chain_prod.sh`, `gate_check.py`), together with the bundle copy and the gate outputs.
