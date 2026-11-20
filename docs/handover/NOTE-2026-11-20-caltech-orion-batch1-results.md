# NOTE 2026-11-20 caltech -> orion (USER): batches 1, 2, 3 + v20g gate, running results

Reply to 111be2f40 / b21808faa / 1997c7d5e / 4d9022943. All four received (35 arms). One running NOTE for all; one line per arm
below as each finishes; per-arm files in `docs/handover/caltech_results/<TAG>/`.

## Setup on our side
- 13 builds (batch 1 x6, batch 2/3 x6, v20g), one worktree per commit, each built on its own H200 node in a holder step.
  Your cmake line, except: **Unix Makefiles instead of Ninja** (no ninja on the Caltech nodes) and a shared ccache launcher.
  Toolchain: cuda/12.2.1 (spack, gcc-11.3.1) + system g++ 11.5 + hpcx 2.17.1 OpenMPI 4.1.7; Kokkos 4.6.02 (d8e9af031) from
  the submodule; `-ffp-contract=off` (the tree's CMakeLists adds `--fmad=false` for device code with it).
- `armc.sh` used as is, with two Caltech changes (copy `armc_cal.sh`): it loads the two modules above (runtime libs), and runs
  `res.py` as a step in the holder (our login-node rule). GPU binding: unchanged (KOKKOS_MAP_DEVICE_ID_BY=mpi_rank; every
  arm is checked for 4 distinct GPUs in mem.log).
- Queue order: P1, P2, P3 + v20g gate (restart lane, hpc-sm-02-16 first), P4, P5, then P6 > P7 (FMP/FMM long pair) > P8 > P9.
  Pairs start together (90 s apart) on comparable nodes, never on 02-16 (suspect for timing; only rt/poison arms go there).
- ETA (rough): builds ~20 min; P1-P3 + v20g within ~1.5 h of that; batch 1 done in ~3 h; batches 2/3 done in ~6-7 h.
  P10 (nsys): will check if nsys is on the nodes; last.

## Results (one line per arm)
