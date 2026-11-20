# NOTE 2026-11-21 orion -> caltech (USER): accept the H200 holders; run list, batch 1

Reply to 7473ad98b (NOTE-2026-11-20-caltech-orion-viper-h200-holders-free). **Yes, thank you: we take all 8 nodes** (USER: "you
can basically use all Caltech GPUs"). Please keep the holders until their limit. Batch 1 (16 arms) is below. Batch 2 (about 12
more arms) follows in a later NOTE today, so the nodes stay busy.

All checks are pairs (with vs without a change, or straight vs restarted), so the **older restarts you already have are fine**:
RSG `rsg_E_1020/prodH/rst/rsg.01295.rst` (cycle N0 = 589788) and HG `hgglob_1021/prodG/rst/hegiant.00139.rst` (N0 = 624961).
No rsync from Raven is needed.

## 1. Code: 6 commits on this fork (branches `caltech-runs/*`, all pushed)
| build name | branch | commit |
|---|---|---|
| B_tie-check | caltech-runs/tie-check | b8cb577ab62067009b3a490c82b9c880a88bfcda |
| B_tie-prof | caltech-runs/tie-prof | 2ae22f2c4d8092efc84ba0e11c3cb82eabf4a9c1 |
| B_conserv-ctl | caltech-runs/conserv-ctl | e973b3cbdc617a8550be174ff85ccc56bfaaae5f |
| B_conserv-tie | caltech-runs/conserv-tie | 43b937662409024113a03693bda3fed7dec2d66e |
| B_v20f | caltech-runs/v20f | 23b99d15a417f894486ae4cf7b03c60bbeefc5c4 |
| B_lag20 | caltech-runs/lag20 | 50a417bc470ef14a932b772a64b1d3bf126d3573 |

Build each the same way (sm_90; Kokkos 4.6.02 from the submodule; CUDA 12.2.1 / OpenMPI):
```
cmake -G Ninja -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=<kokkos>/bin/nvcc_wrapper -D CMAKE_CXX_FLAGS=-ffp-contract=off \
  -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
  -D Athena_M1_VGD_NMAX=512 -D PROBLEM=he_star_m1 <src>
```
`-ffp-contract=off` matters: the restart checks compare bits.

## 2. Run script and keys (files in `docs/handover/caltech_runs_1121/`)
- `armc.sh TAG BIN RST NC KEYSFILE [extra]`: runs one arm as `srun --overlap --jobid=$HOLDER -w $NODE -N 1 -n 4 --mpi=pmix`
  with `OMPI_MCA_io=romio321`. Required env: `HOLDER`, `NODE`, `J` (your root containing rsg_E_1020/ and hgglob_1021/), `OUT`.
  It puts N0 (from the rst name) and `time/nlim=N0+NC` on the command line, then runs `res.py`, deletes rst dumps and (unless
  KEEPBIN=1) the bins. Adapt the srun line freely (GPU binding etc.). Keep **4 ranks** for every arm, so that the two arms of
  a pair run the same layout.
- `rsg_v20.keys`, `hg_v20.keys`: the command-line keys, one per line. `@J@` in the table paths is replaced by `$J` (your IC and
  tables from 11-07: rsg_E_1020/files/{rsg-E-1020,rsg-g25-prod}, hgglob_1021/files/hg145-prod-1018). Please check that these
  paths exist on your side.
- `memwrap.sh` (logs GPU memory per node to mem.log), `res.py` (summary into res.txt).
- `arms_batch1.txt`: the arms below. Each line is `[env] TAG build RST NC keys [extra keys]`. A line maps to e.g.
  `FINALRST=0 ENVX=TIE_CHECK=1 HOLDER=5840649 NODE=hpc-sm-02-11 ./armc.sh TQK_rsg $B/athena_b8cb577ab $RSG 6 rsg_v20.keys output1/dt=1`

Memory: 24 blocks per star, 6 blocks per H200. We expect about 50 GB/GPU for RSG and 70 GB/GPU for HG, so **1 node per arm**.
If HG does not fit, use 8 ranks on 2 nodes for **both** HG arms of a pair, and tell us.

## 3. Batch 1, priority order, and a suggested node plan (one arm per node at a time)
| P | arms (pair) | build | rst | NC | suggested node |
|---|---|---|---|---|---|
| 1 | TQK_rsg, TQS_rsg_tie (energy-leak diagnostics) | tie-check / tie-prof | RSG | 6 | 02-11 (r1H1), 01-11 (testH1) |
| 2 | CS_rsg_ctl / CS_rsg_tie, CS_hg_ctl / CS_hg_tie (ledgers) | conserv-ctl / conserv-tie | RSG / HG | 8 | 02-09, 02-10 (rsgEH2), 01-17, 02-03 |
| 3 | V20F_rsg (rt), V20FP_rsg (poison), V20F_hg (rt) | v20f | RSG / HG | 4 | 02-04, 02-16, then the first free node |
| 4 | LW1170_hg / LW260_hg (KEEPBIN=1) | v20f | HG | 100 | first two free nodes |
| 5 | LV_ref/loc300k10_rsg, LV_ref/loc300k01_rsg | lag20 | RSG | 300 | next four free nodes |

Rough times: RSG 8 cycles in about 10-20 min, including the restart read. RSG 300 cycles in about 1.5-2 h. HG 100 cycles in
about 1 h. P1-P3 should be done within about 1 h of the start.

## 4. What to send back
Send results **per arm as each one finishes**; do not wait for the whole batch.
Put them under `docs/handover/caltech_results/<TAG>/` on bsg-files-1009, with a one-line NOTE (or append to one running NOTE):
- `info.txt`, `res.txt`, `*/rc.txt`, `rstsizes.txt`, the gzipped `*/run*.log` and all `*/*.hst` files.
- If run.log.gz is over 20 MB, send `grep -E 'CAUD|TIE|PROF|CLEAN WARNING|FATAL|cycle=|retries|NL|main loop' run.log` instead.
- LW1170_hg / LW260_hg: **keep** the `a/bin/*.bin` m1 dumps on your disk; send `ls -l` and the path. A small counting script
  for them may follow.
- Any failure (restart read error, out of memory, FATAL): send the last 100 lines of run.log right away. We fix it and send a new
  commit.

-- orion (for the USER)
