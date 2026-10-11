# TASK raven -> DeltaAI: BSG true repro with the recommended scheme, R vs RC smoke + cycle-0 IC check (1 short ghx4-interactive job)

From: Raven (worker D, bsgrec-1011). The Raven gpudev slot is ~1 h out and only 1 node (4 A100 x 40 GB: 4 blocks/GPU do not fit),
the Raven gpu partition is ~1 day. **First result wins**; Raven cancels its gpudev smoke 31078483 when you post "STARTED".
Short test only: <= 1 h wall, 1 node. Production of these arms runs at Caltech (TASK-2026-10-11-caltech-bsg-rec-consrc.md);
this smoke tells Caltech whether the decks are clean and how much memory hesdirk2 needs.
Same environment rules as your earlier BSG tasks: `MPICH_GPU_SUPPORT_ENABLED=1`, plain srun, 1 rank per GPU.

## 1. Build (report md5)
- Fork branch **`prodfinal-1011`** @ **530c38e2** (if the tip has moved, build the tip and report the hash; the decks name every
  key whose default changed, so the result is the same). `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1, GH200), as your
  6c5d8fb2 / d28dacd8 BSG builds.

## 2. Files
```
git fetch fork bsg-files-1009
git archive fork/bsg-files-1009 docs/handover/bsg-hrdet-1009 docs/handover/bsgrec-1011 | tar -x -C $G/bundle
bash $G/bundle/docs/handover/bsg-hrdet-1009/SETUP.sh $G/files      # or reuse your 10-09 files dir (md5 check)
for a in R RC; do sed "s#@FILES@#$G/files#g" $G/bundle/docs/handover/bsgrec-1011/bsg$a.athinput.in > $G/bsg$a.athinput; done
```
- `bsgR` = the 10-09 BSG production input `bsg3d_truerepro2_hr_lm.athinput` (full 256^3 mesh, 16 MeshBlocks 256x64x64,
  `vet_gd_twin_lowmem = true`) with the recommended scheme (hesdirk2, time2_vet_col predict, vimp, vstage, xthin 0, fse, hc), the
  current defaults named (vet_seg emit, face, ...), `implicit_gas_src = matched`, and the ORIGINAL IC named:
  `ic_rad = none`, `he_ic_rad = file`, `he_ic_rad_gas = column`, `he_ic_rad_hse = false`.
- `bsgRC` = bsgR with ONLY `implicit_gas_src = conservative`.
- Command-line keys must already be in the deck (they are: `time/nlim`, `time/ndiag`).

## 3. Runs (one job, 1 node x 4 GH200, ghx4-interactive, <= 1 h; stop the job if anything below FATALs)
| run | ranks x blocks/GPU | input | extra keys |
|---|---|---|---|
| R4  | 4 x 4 | bsgR  | `time/nlim=10 time/ndiag=1` |
| RC4 | 4 x 4 | bsgRC | `time/nlim=10 time/ndiag=1` |
| R2  | 2 x 8 (Caltech layout, memory only) | bsgR | `time/nlim=4 time/ndiag=1` |

Each in its own run dir (`-d`). The cycle-0 dumps (bin/bsg3d.hydro_w.00000.bin, bsg3d.m1.00000.bin) are written at start; keep them
for R4 and RC4. If R2 runs out of memory, just report that (the 10-09 be/hr peak was 65.0 GiB at 8 blocks/GPU, 42.9 GB at 4).

## 4. Report: NOTE-2026-10-11-deltaai-bsgrec-smoke.md (push to bsg-files-1009)
- **STARTED** line as soon as the job runs (job id), so Raven cancels its copy.
- Per run: rc, FATAL / NaN / NON-CONVERGED counts, GPU peak per GPU (nvidia-smi MiB), s/cycle (steady, cycles 3-10, and the
  rotation cycles separately), dt, Picard mean/max (run.log `implicit transport: solves=... Picard iterations mean=... max=...`),
  any `NOT ADMISSIBLE` / stage fallback / `time2_pred_cap ... limited` lines, the `implicit_gas_src` banner if printed, setup time.
- From user.hst of R4 and RC4: L_top, L_in, Etot per cycle (paste the 11 rows of t, L_top, L_in, Etot, Picard for both).
- **Cycle-0 IC check** (the user wants the ORIGINAL 10-09 IC): for R4 and RC4 run
  ```
  PYTHONPATH=$G/bundle/docs/handover/bsg-hrdet-1009/page python3 $G/bundle/docs/handover/bsgrec-1011/ic0check.py \
      check <run dir> $G/bundle/docs/handover/bsgrec-1011/ic0_ref_1009.json
  ```
  `ic0_ref_1009.json` is the fingerprint (per variable: md5 of the assembled float32 array, sum/min/max, shell means) of the
  Raven 10-09 hr cycle-0 dumps (binary 09872817 = hrup a5de361a + truerepro2 d385de40; the d385de40 cycle-0 dumps are bitwise
  identical to it, all 9 variables). Paste its output: BITWISE per variable, or the max relative differences (GH200 vs A100 device
  round-off would be ~1e-7; anything larger means the IC is not the original).
- GPU-h used.
