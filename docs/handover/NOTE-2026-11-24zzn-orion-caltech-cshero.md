# NOTE 2026-11-24zzn orion -> Caltech: CSHERO cubed-sphere imprint arms A / B / C on bsg's node when bsg reaches tlim (~16:45 PST)

Run these only on GPUs freed by bsg's end (or any GPUs not used by production). Never slow the production runs. Files are in docs/handover/cshero/.


**Why.** The docs hero image (cs RSG, cycle 1200, t = 13.89 d, rho/<rho> - 1) shows a 4-fold X along
the cube's panel diagonals (~3e-4) and a rippled ±1.5e-3 ring at the photosphere. It came from the
RELEASE template (38f3d3c95). That template already names cs_vertex_mono = true, which is also the cs
default. It ran the K = 4 window ring (vet_rotate_avg_ema = 0), which has no well-balanced form.
T2 (262affd6f / 41a7b4952) switched the template to the EMA ring (vet_rotate_avg_ema = 128 +
vet_rotate_avg_ema_wb = true). These three arms test whether that switch removes the X.

**Binary.** athena_ia123rc4v_262affd6f (md5 7d7c4a71), the he_star_vet build. The template needs
`-D PROBLEM=he_star_vet`, so this binary is the right one.
- Between 262affd6f and 41a7b4952, `git diff` over src/ and over the template is empty (only tests/docs
  changed).
- Please write the binary md5 into each run dir.

**Decks** (copy these into docs/handover/cshero/). All are the 262affd6f template
inputs/rad_vet_clean/star3d_cs.athinput, with nlim = 1200 and hydro_w + vet dumps every 150 cycles:
| file | arm | change vs the template |
|---|---|---|
| star3d_cs_A_hero.athinput | A: the hero's closure | vet_rotate_avg_ema = 0, no vet_rotate_avg_ema_wb (K4 window) |
| star3d_cs_B_prod.athinput | B: production keys (T2 template as is) | none |
| star3d_cs_C_fixed.athinput | C: fixed-set control | B with vet_gd_rotate_dt = 0, vet_gd_rotate_every = 0 |

**Run.** Use one dir per arm, and link star3d_data from the 262affd6f checkout
(inputs/rad_vet_clean/star3d_data; the deck reads star3d_data/rsg_*.txt).
- 4 GPUs per arm, as for the hero (24 MeshBlocks of 548x16x16).
- On an 8-GPU node, run A + B together, then C. Or run all three at 2 GPUs each.
```
ln -s $ATHENAK/inputs/rad_vet_clean/star3d_data star3d_data
mpirun -np 4 -x KOKKOS_MAP_DEVICE_ID_BY=mpi_rank -x OMP_NUM_THREADS=1 $EXE \
  -i star3d_cs_<ARM>.athinput -d run_cs > run_cs.log 2>&1
```
- Fresh start only. Do not add keys on the command line: a key named only on the command line of a fresh
  start is FATAL for vet_rotate_avg_ema_wb.
- No restart dumps are written.

**Analysis.** cshero_imprint.py (numpy only):
```
python3 cshero_imprint.py run_cs/bin/rsg.hydro_w.*.bin > imprint_<ARM>.txt
```
- It prints a table per dump: per radius bin, the maximum of:
  - the equatorial m = 4 and m = 8 Fourier amplitudes of rho/<rho> - 1 (the X at the seams = m = 4);
  - sym_rms: the panel-symmetric part of drho, i.e. drho averaged over the 6 panels at each lateral
    index, minus its solid-angle mean (grid imprint; convection averages out);
  - rms_all: the total rms.
- It also writes `<dump>.imprint.npz` next to each dump.

**Post** (NOTE + push):
- imprint_A.txt, imprint_B.txt, imprint_C.txt and the .npz files (small);
- for each arm: run_cs.log tail (s/cycle, CLEAN WARNING / FATAL counts), user.hst, and the binary md5;
- if possible, the hydro_w dumps at cycles 0, 600 and 1200 of B, so the hero can be regenerated from B
  (6 x ~160 MB per arm at double precision; skip if too large, and tell us);
- delete the other dumps after the analysis.

-- orion (coordinator)
