# TASK for DeltaAI: He giant fine-grid run (N897), smoke + production chain (duplicate of viper hegiant897)

## 1. What / why

- **What:** 3-D radiation-hydro (implicit M1, vet_gd closure, TOPS opacities, general EOS table) of a 2.56 Msun He
  giant envelope before core collapse (MESA profile), spherical wedge r 3-200 Rsun x (pi/3 x pi/3), 897 x 128 x 128.
- **State:** restart `hegiant.00022.rst` = the viper scout128 run at t = 9.004208e5 s (10.42 d, cycle 47861), remapped
  onto the N897 radial grid (fine band 60-95 Rsun); TLIM 2.592e6 s (30 d).
- **Why here:** the viper chain `hegiant897` (jobs 12119930-39, 2 nodes x 2 MI300A, 6 h links) is still PENDING.
  Duplicates are intended: **leave the viper copy alone, the user decides** which copy continues. Do not touch
  any other run.

## 2. Code

- Branch **`hegiant-opn-1007`** on the fork, commit **f3a66907** (`git fetch fork hegiant-opn-1007`; the commit adding
  this note sits on top and changes docs only). f3a66907 = hegiant-1006 (he_star_m1 pgen: he_base_heat,
  he_gm_column, MLT scaffold ramp) + opacnewt-cliff-1007 (`rad_m1/implicit_opac_newton_slope_max`). It contains the
  CUDA fixes 952abbfe (c2p_track lambdas) and d3c37ba6 (CoordData init).
- viper binary: `athena_he_gpu72_f3a66907_hegopn`, md5 fc7b33a1ea517e9bddd83850f5b2d0f6 (ROCm 7.2, MI300A).
- Build on DeltaAI with `build_inc_deltaai.sh`, new target PROBLEM=**he_star_m1** (GPU; a CPU binary is not needed
  for this task). Stack as the BSG task: PrgEnv-gnu, gcc-native/14, cudatoolkit/25.5_12.9, cray-mpich/9.0.1,
  craype-accel-nvidia90, `Kokkos_ARCH_HOPPER90=On` + `Kokkos_ARCH_ARMV9_GRACE=On`, nvcc_wrapper, MPI on. Record md5.
- **FMA:** nvcc contracts a*b+c into FMA by default, the MI300A build differs in where it does
  (NOTE-2026-10-02-deltaai-coorddata-init.md: a Picard pass-count flip was FMA alone; Grace `-ffp-contract=off`
  reproduced x86 exactly). So expect agreement with viper to ~1e-6 relative over 60 cycles, not bitwise, and
  occasional +-1 Picard passes. Do not build with contraction off for production (speed); use it only if a smoke
  disagrees beyond the tolerances of section 5 and you need to separate FMA from a port bug.

## 3. Bundle

The user copies `/viper/ptmp2/jinma/hegiant_deltaai_bundle` (2.2 GB, almost all the restart) to DeltaAI, e.g.
`scp -r viper:/viper/ptmp2/jinma/hegiant_deltaai_bundle <deltaai>:/work/nvme/bivj/jma20/` or Globus. Then:

```
bash /work/nvme/bivj/jma20/hegiant_deltaai_bundle/SETUP.sh /work/nvme/bivj/jma20/hegiant_deltaai_bundle
```

SETUP.sh checks all 21 files against `MD5SUMS` (restart md5 cdae07fcaef0d942df8d1148c730dbcf) and writes `local/`
with the bundle path filled in. **Keep the bundle where it is for the whole run**: the restart embeds the viper paths
of three files that he_star_m1 re-reads on every restart, and `local/restart_keys.txt` overrides them:

| file | what | md5 |
|---|---|---|
| `rst/hegiant.00022.rst` | remapped N897 restart, t 9.004208e5 s, cycle 47861, 4 MeshBlocks | cdae07fc... |
| `ic/ic_giant_own.txt` | IC column (`problem/he_ic_file`: HSE / a_ref / frozen MLT scaffold targets) | ddc72692... |
| `tables/rosseland_tops_hegiant_blend.txt`, `planck_tops_...` | `problem/he_opac_table`, `he_planck_table` | c1ef063c..., b699c60d... |
| `local/hegiant_n897_rst_embedded.athinput` | the input embedded in the restart (reference; the restart carries it) | |
| `local/hegiant_scout128_N445_fresh.athinput` | viper run.cfg `IN` (N445 fresh start; not used on restart) | |
| `local/dual_hegiant_deltaai.sh` | chain link script (port of the viper one, see section 6) | |
| `local/check_link.py`, `local/ana_run.py` | run-dir summary: FATAL/NON-CONV, s/cycle, dt, Picard, R_ph, max v_r | |
| `scripts/rst_info.py`, `scripts/linkcheck.py` | restart time/cycle + output last_time keys; pre-link sanity check | |
| `smoke_ref_viper/` | viper smoke 12120235: run.log, hst files, job output, run.cfg | |

Verified on viper: the restart-embedded input names no other file (grep of the header up to `<par_end>`: only the
three problem/ paths; EOS table is built in memory, grid is the analytic r-stretch polynomial; `he_grid_dump` writes
`hegiant.x1grid.txt` into the run dir). `check_link.py` needs numpy, scipy, h5py.

## 4. Layout

4 MeshBlocks of 897 x 64 x 64 (mesh 897 x 128 x 128, 2 x 2 blocks in theta/phi). DeltaAI: **1 node x 4 GH200, 4 ranks,
1 block per GPU** = the viper production layout (2 nodes x 2 MI300A). The viper smoke used 1 node x 2 GPUs (2 blocks
per GPU), so its s/cycle is NOT the per-GPU production load: expect roughly half per cycle with 4 GPUs.
Launch recipe (as the BSG task): partition ghx4, every rank sees all GPUs, `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank`,
`MPICH_GPU_SUPPORT_ENABLED=1`, `OMP_NUM_THREADS=1`, `srun -n 4 -c 16 --cpu-bind=cores`; NOT
`--gpus-per-task=1 --gpu-bind=closest`.

Restart command (link 1; the link script builds it from run.cfg):

```
srun -n 4 -c 16 --cpu-bind=cores <BIN> -r <run>/rst/hegiant.00022.rst \
  output1/last_time=900400 output2/last_time=885600 output3/last_time=885600 output4/last_time=900000 \
  output5/last_time=864000 output6/last_time=885600 time/tlim=2.592e6 \
  time/cfl_number=0.3 time/restart_refill_ghosts=true rad_m1/implicit_opac_newton_slope_max=3 \
  $(cat <bundle>/local/restart_keys.txt) -d <run> -t <wall - 20 min>
```

(`output*/last_time` come from `python3 scripts/rst_info.py <rst>` line 2; they stop catch-up dumps.)

## 5. Smoke (60 cycles, before any production; again after ANY change of binary, keys or run.cfg)

```
mkdir -p <smoke>/rst && ln -s <bundle>/rst/hegiant.00022.rst <smoke>/rst/
cp <bundle>/local/run.cfg.template <smoke>/run.cfg      # fill BIN, BIN_MD5, ENVSH
sbatch -J hegiant_smk -t 00:20:00 -o <smoke>/dual.%j.out -e <smoke>/dual.%j.err \
  <bundle>/local/dual_hegiant_deltaai.sh <smoke> 60
python3 <bundle>/local/check_link.py <smoke>
```

The run never writes into the bundle (new rst go next to the symlink in `<smoke>/rst`; the smoke ends with a 4.25 GB
`hegiant.00023.rst` written at nlim).

Viper reference, smoke job **12120235** (same restart, same keys, binary md5 fc7b33a1, 1 node x 2 MI300A,
`he_giant_1006/runs/smk_smax3`, copied to `smoke_ref_viper/`):

| quantity | viper |
|---|---|
| rc / FATAL / NaN | 0 / 0 / 0 |
| dt at cycle 47861 (first step after restart) | 9.223353e+00 |
| dt at 47862 / 47871 / 47881 / 47921 (last) | 1.844671e+01 / 1.856086e+01 / 1.856955e+01 / 1.858862e+01 |
| Picard (end-of-run `implicit transport:` line) | solves 60, mean 18.5, max 42, NON-CONVERGED 0 |
| `NEWTON-FALLBACK` log lines (60 cycles) | 960 (informational; scattered cells 4.8-11 e12 cm) |
| last hst row (t, dt, mass, tot-E) | 9.0152570477933483e+05, 1.8588619517698017e+01, 1.5435375539788510e+32, 1.5044198258655804e+47 |
| `he_ic_balance` line at startup | max \|rho/rho_col - 1\| = 0.00515103 at r = 2.15997e+11 (proves the IC file was read) |
| check_link.py | R_ph(tau_R 2/3) 73.316 Rsun, r(rho 1e-12) 76.048 Rsun, max \|<v_r>\| 37.579 km/s at 72.79 Rsun (t 10.434 d) |
| s/cycle (2 blocks per GPU) | 2.779 (check_link.py, cycles 6-60) |

**Pass:** rc 0, 0 FATAL/NaN, 0 NON-CONVERGED; the same `he_ic_balance` numbers (to the printed digits: proves the
bundle files are read); dt at the listed cycles within 1e-5 relative; last-row mass and tot-E within 1e-6
relative; Picard mean within ~10 % (16.7-20.4); R_ph to 0.01 Rsun. Report s/cycle (median of cycles 10-60 from the
`elapsed=` lines). A disagreement at 1e-3 or worse in mass/energy or dt is a port bug: stop and report, no production.

## 6. Production chain

- Run dir `<run>` (under `/work/nvme/bivj/jma20/`), `rst/` holding a symlink to the bundle restart, `run.cfg` from the
  template. Link script `local/dual_hegiant_deltaai.sh <run>` (no 2nd argument): guard against an older running job
  of the same name, CANCEL/DONE/STOP files, binary md5 check, restart from the newest `<run>/rst/*.rst` with the
  last_time keys, linkcheck of the previous run.log (NaN/FATAL/error, a newer rst, dt collapse < 1e-3 x median),
  code wall = job limit - 20 min, STOP on FATAL/NaN or rc != 0, DONE at TLIM. Fill the `#SBATCH -A` line.
- TLIM **2.592e6 s** (30 d); wall limit per link = the ghx4 maximum (`sinfo -p ghx4 -o %l`), links chained with
  `-d afterany:<prev>` (the script exits 0 without running when it should not). Remaining 1.69e6 s at dt ~18.6 s =
  ~91,000 cycles; the cost per cycle on 4 GH200 comes from your link-1 log (not measured on viper with 4 ranks).
- Output cadence (from the restart, unchanged): hst every 25 s, log 1e4 s, bin hydro_w / m1 / m1_vet every 0.25 d
  (~765 MB per set), rst every 0.5 d (~4.25 GB each): ~80 bin sets + ~40 rst to 30 d (~230 GB); thin old rst if the
  quota needs it, keep the newest two and every 5 d.
- **Drop `time/restart_refill_ghosts=true` after the first link** (it is only for the remapped restart): edit XKEYS in
  `<run>/run.cfg` while link 2 is still pending, run a 10-cycle smoke of the newest rst with the new run.cfg in a
  scratch dir (`dual_hegiant_deltaai.sh <scratch> 10`), and keep link 2 queued only if it passes (else hold it).
- Do not change any other key mid-run. The physics schedule (MLT scaffold ramp ended at 10 d) is in the restart.

## 7. Report back

`docs/handover/NOTE-2026-10-0x-deltaai-hegiant.md` on branch `hegiant-opn-1007` of the fork (fetch, merge, push;
never force): binary sha + md5 + stack; smoke job id and the section 5 table side by side with viper; production
job ids per link, s/cycle, dt, Picard mean/max, NON-CONV count, simulated time per link, run dir. The full hst files
grow to ~60 MB at 30 d (852 bytes per 25 s row): push them thinned to every 40th row plus the header
(`awk '/^#/ || NR % 40 == 0'`, ~1.5 MB) and the `check_link.py` output with the NOTE when a link ends; no interim
reports. **Do not touch other runs; the user decides about the viper copy (hegiant897).**
