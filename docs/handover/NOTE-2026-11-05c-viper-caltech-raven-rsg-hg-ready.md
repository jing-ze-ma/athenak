# NOTE 2026-11-05c viper -> Caltech: RSG E + HG on Raven, items 1-3 READY (reply to TASK 7b602365e)

Nothing has been submitted. A launch waits for: the rst copied by the user, then a passing smoke, then SMOKE_OK.

## 1. Build
- Binary: `/raven/ptmp/jinma/raven_gap_1105/bin/athena_he_a100_nofma_23a55658_fmadoff`, **md5 620c85529afb5be4033203701c59b0d0**.
- Source: exactly 23a5565803c8 (de528da5 is an ancestor, so no patch). Kokkos d8e9af03, as Caltech.
- Recipe = Caltech's:
  - PROBLEM=he_star_m1, MPI on, CUDA, Release;
  - `-ffp-contract=off --fmad=false`;
  - AMPERE80 instead of HOPPER90; gcc/13, cuda/12.6, openmpi_gpu/5.0.
- The only other difference: `Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF`, as in all our A100 he builds. It is an allocator setting with no effect on results.
- `athena -c` OK on gpudev.

## 2. Space and copy destinations
- /raven/ptmp has 470 TB free (97 % full). The per-user quota is unreadable: mmlsquota reports "GPFS down" on the login nodes, and no quota/lfs tools exist.
- Rough output estimate for 4 days: HG ~0.9 TB, RSG ~2.5-4 TB. rst files are pruned to the newest 3 per run.
- **Copy ONE rst per run** to:
  - RSG E: `/raven/ptmp/jinma/raven_gap_1105/rsgE/rst/` (rsg.NNNNN.rst);
  - HG: `/raven/ptmp/jinma/raven_gap_1105/hg/rst/` (hegiant.NNNNN.rst).
- Also copy the Caltech hst files to `.../{rsgE,hg}/caltech_hst/` (for the smoke comparison), and the rst md5s.

## 3. Run dirs and launchers
- `run.cfg` KEYS are the Caltech last command lines verbatim, without -r/-d/-t and last_time. viper checked this programmatically for both runs.
- The last_time values are recomputed from the rst that arrives, as floor(t/dt)*dt, honouring the outputN/dt overrides.
- `jobs/gaplink.sh`:
  - layout: 6 nodes x 4 A100 = 24 ranks x 1 block;
  - restarts from the newest rst, `-t 23:30:00` per 24 h link;
  - checks the binary md5;
  - **STOP** on FATAL / DIVERGED / final NON-CONVERGED > 0 / nan / rc != 0;
  - CANCEL / DONE / FINAL control files.
- `jobs/submit_chain.sh <rsgE|hg> [NLINKS]` chains the links with afterany and refuses to run without SMOKE_OK.

## 4. Smoke plan
- `jobs/smoke.sh <rsgE|hg>`: 30 cycles on 6 nodes.
- Pass:
  - rc 0, no FATAL/nan, NON-CONVERGED 0;
  - GPU memory < 40 GB/A100 (Caltech ~19 GB/block for RSG, ~25 GB/block for HG);
  - `jobs/compare_hst.py`: the first hst row agrees with the Caltech row at the same time to ~1e-12..1e-10 relative. Above 1e-8: no start.

## 5. Reverse plan
- When the Caltech holders start, `touch FINAL` in the run dir. The running link then ends at its wall-limit rst and the chain stops.
- For an earlier stop, shorten the next link's TimeLimit with scontrol.
- There is no signal handler: a plain scancel loses everything back to the last rst.
- The user copies the newest rst back to prodH/rst (RSG) or prodG/rst (HG).
- Details: `/raven/ptmp/jinma/raven_gap_1105/README.md`.
