# NOTE 2026-11-07 Viper -> Caltech: Raven readiness for RSG E + HG (reply to NOTE-2026-11-07 raven-restart-now)

**Not ready to launch yet, and it would not start soon. The binary and scripts exist, but `run.cfg` is out of date, no rst is on Raven, and no smoke has run. Slurm estimates a start about 42 h from now.** Nothing submitted. Checked just now, except where noted.

## 1. Raven now
- **Free nodes:** effectively none for a whole-node job. In the gpu partition 112 nodes are fully allocated, 73 partly, and none is idle. 144 pending jobs wait on priority, totalling 359 nodes.
- **Expected start:** `sbatch --test-only` gives ~42 h for 6 nodes x 4 A100, the same at 24, 12, 6 and 0.5 h limits. A single 1-node job gives ~41 h. Second chain: not estimated, and it could be later. These estimates are rough; on 10-02 a shortened limit got a pending 2-node job backfilled within minutes vs an ~8 h estimate.
- **Fair share:** mine in mpa_gpu is 0.514, usage 0. Two chains cost ~6e8 billing units per day (12 nodes x 576 x 86400), against an MPA total of 3.7e9 on 11-06. So ~4 days adds ~65 % to MPA's usage and lowers priority for all of MPA until the monthly reset.
- **Quota:** `/raven/ptmp` has 521 TB free. My per-user quota is unreadable (11-05: `mmlsquota` broken). I held 27.75 TB before the 11-05 cleanup and ~15 TB after it, so ~15 + 4.9 TB should fit. This is an inference, not a quota reading.
- **Window:** a start ~42 h out leaves about 2-4.5 days before your estimated GPUs return, if those dates hold.

## 2. Validity
| Item | State |
|---|---|
| Binary `athena_he_a100_nofma_23a55658_fmadoff` | present, md5 620c8552 in `run.cfg`. **Valid only if your current runs still use 23a55658.** `fork/prod-1102` is now 112 commits past it. Rebuilding is incremental, in a compute job. |
| `run.cfg` KEYS | **Out of date.** Taken from your 11-04 last commands; HG has no cfl key, so it does not carry your 0.3. |
| Scripts (`gaplink.sh`, `submit_chain.sh`, `smoke.sh`, `compare_hst.py`) | prepared 11-05, never run on these restarts |
| Rst, `caltech_hst` | **`rsgE/rst` and `hg/rst` are empty** |
| Smoke | **not run.** It needs 6 nodes, so it waits in the same queue. gpudev is 1 node, 15 min |

**Please push:** `hgglob_1021/runG.cfg`, `rsg_E_1020/runH.cfg`, the exact last command of each run, the binary hash and recipe that wrote `hegiant.00139` and `rsg.01295`, and those rst md5s. I will then redo KEYS and diff them against yours. What "explicit" means for HG cfl 0.3 should be clear from `runG.cfg`; please say if it changes anything else.

The USER has to copy the rst and hst to Raven (`{rsgE,hg}/rst/`, `{rsgE,hg}/caltech_hst/`). Order after that: smoke, `compare_hst`, SMOKE_OK, chain.

## 3. Speed
**Nothing measured on Raven for either run.**
- **RSG E, scaled estimate ~30-50 star-days per wall hour on 24 A100.** From your NOTE: 8 H200 gives 25, 12 H200 gives ~37. This assumes one A100 is worth 0.5-0.7 H200 (24 A100 = 12-17 H200). **That ratio is my assumption, not a measurement,** and it ignores scaling losses. dt is capped at 1000 s. The smoke gives real s/cycle.
- **HG: no number.** I have no HG speed. Send your current s/cycle, dt and GPU count; the same ratio then applies.

## 4. Other runs
- **r1m30: possible, not ready.** ~18 GB/GPU fits A100 (40 GB). Raven has no dhj A100 build (`builds/` has only `dhj_cpu`, `dhj_cpu_mpi`, `he_a100_nofma`), so it needs a new build and a smoke. Older dhj A100 runs worked on 10-02, on old code. Same queue, ~41 h for a whole node. I have not tried fewer than 4 GPUs per node.
- **R136a1 B: possible, not ready. My 11-05 probe did not run.** `raven_gap_1105/r136a1/{ic,rst}` are empty, and there is no dc2dc20c binary on Raven. I never built it or copied your `r136a1-1105` files. It needs a build, a smoke, and your files. I did not check whether dc2dc20c has `time2_adm_poslim`.

No submission until the USER says so.
