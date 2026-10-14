# TASK raven -> DeltaAI: VALIDATION ONLY of fin4-linhao-bc (Linhao IC + Linhao's exact inner BC) on 1 node x 4 GH200

From the Raven coordinator (worker DA). USER: DeltaAI is for short tests only. This is NOT a production run: no chain, no GO file,
no links after the gate. It is an independent GPU validation for worker XB (linhaobc_1014), whose Raven/Viper dev-node checks are queued.
One job, <= 45 min wall (ghx4-interactive or ghx4, whichever starts sooner; say which and the estimate in the NOTE).

## 1. Build (report md5)
- Fork branch **`fin4-linhao-bc-eb1e032c` @ `eb1e032c`** (= fin4-linhao d32aebaf + 70233620 MB ghost bath + XB's opt-in `problem/he_bc_linhao`).
  It is a snapshot of XB's branch; XB is still working, so a later NOTE may name a newer hash. Build exactly eb1e032c for now.
- `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1), same flags as your 9ceb9cdc build. Name `athena_hes_gpu_eb1e032c`.
- Report the md5 and `git rev-parse HEAD`.

## 2. Files (`docs/handover/bsgli-val-1014/`, check MD5SUMS)
- One dir `<F>` with:
  - your 10-09 BSG files (ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805; the deck still reads them);
  - `linhao_profile_full.txt` = gunzip of `linhao_profile_full.txt.gz` (17 MB gz; the unpacked md5 must be 41edc225).
- Deck `bsgLI_bc_deltaai.athinput.in` -> `sed "s#@FILES@#<F>#g"`. Change no other key.
  It is Raven's `prod/LI/bsgLI_full_fin3_bc.athinput` (md5 d37e861a) with only three changes:
  - the 4 file paths are templated;
  - `vet_gd_halo_ipc = true` and `vet_gd_twin_lowmem = true`, your single-node BSG layout keys. Both are bitwise, and lowmem takes ~50 instead of ~79 GB/GPU.
  The key under test is `<problem> he_bc_linhao = true`.
- Scripts: `val.sh` (the job), `rst_pick.py`, `chaincmp.py`, `cmplast.py`, `bin_convert.py` (the e1caltech-1013 copies).

## 3. Layout
1 node x 4 GH200, 4 MPI ranks x 4 MeshBlocks (256 x 64 x 64), 1 rank per GPU, OMP_NUM_THREADS=1, MPICH_GPU_SUPPORT_ENABLED=1.
This is your bsgprod/bsgaccel N4 layout. Use your usual launch (`LAUNCH=` in val.sh, e.g. with your gpuwrap.sh).

## 4. Run (one job, scratch dir, nothing else touched)
```
sbatch --export=ALL,W=<abs dir of bsgli-val-1014>,BIN=<bin>,DECK=<resolved deck>,S=<scratch dir>[,LAUNCH="..."] <your flags> val.sh
```
- **C**: fresh, 60 cycles.
- **L1 + L2**: L1 runs fresh to cycle 30 and ends with the link-end rst. L2 is a NEW srun from rst_pick.py's pick, with its last_time keys, to cycle 60.
- **Coff**: C with `rad_m1/vet_gd_halo_ipc=false`.
- Test-only keys on the command line: ndiag 1, hst every 200 s, hydro_w/m1 bins every 1000 s, the other outputs off.
- **Gate**:
  - `chaincmp.py C L1 L2` -> CHAINCMP BITWISE: every hst row and every scheduled dump present once, data arrays bitwise;
  - `cmplast.py C L2` and `cmplast.py C Coff` -> CMPLAST BITWISE;
  - Coff hst rows BITWISE.
- If the gate is not bitwise, do not rerun with changed keys: report.

## 5. Report: NOTE-2026-10-14-deltaai-bsgli-val.md (push to bsg-files-1009)
- Binary md5 + HEAD, the deck md5 after the sed, the profile md5, job id, partition and the queue estimate.
- Per arm: rc, FATAL / NaN / DIVERGED counts, the last `Picard iterations mean= max= NON-CONVERGED=` line, the dt values, and the halo_ipc line.
- The he_bc_linhao and he_ic_linhao log lines.
- **s/cycle** over cycles 10-60 of C, the init time, and the **GPU memory** per GPU (peak).
- The gate output: CHAINCMP / CMPLAST lines, the Coff verdict, and the number of hst rows and dumps compared.
- Keep the scratch dirs. Push C's `bsg3d.user.hst` + `bsg3d.hydro.hst` (small) so XB can compare them with Raven A100 and Viper MI300A.
