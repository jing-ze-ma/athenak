---
name: exo-fms-ck-tables
description: Where the Exo-FMS correlated-k tables live now, and the two traps hit while installing them (GPFS quota truncation, ROMIO on /tmp)
metadata:
  type: project
---

Fetched 2026-08-23 from `https://github.com/ELeeAstro/Exo-FMS_column_ck` (`git clone
--depth 1`, 92 MB; only ~3.1 MB of it is needed). Files verified against the reader in
`deep_hot_jupiter_rt.cpp`: `38 34 11 8`, T = 100..6100 K, p = 1e-8..1000 bar, band edges
DESCENDING 324.68 -> 0.26 um, exactly as `PROVENANCE.md` claims.

## Where they are

`data/exo_fms_ck/` **in the repo** — the single copy, installed 2026-08-23 after the
quota cleanup and verified by md5 on all 15 files. Gitignored, so it is not tracked.
Staging copies on `/freya/ptmp` and `/tmp` were deleted once this one was verified.

The first attempt at this copy hit a full quota and left **zero-byte stubs** rather than
failing loudly (`cp` printed the error but kept the truncated files) — the same
silent-truncation trap as [[gpfs-quota-wall]]. Always `find data/exo_fms_ck -type f -empty`
after any re-copy.

Run with `problem/ck_data_dir` and `problem/ck_table` pointing at absolute paths under
`data/exo_fms_ck` — the doc notes they resolve against the WORKING directory, which `-d`
changes. Output goes to [[test-output-location]], never next to the tables.

## ROMIO cannot write AthenaK output on /tmp

A run with `-d /tmp/...` dies with a bare `MPI_ERR_NO_SUCH_FILE: no such file or
directory` and `MPI_ABORT`, with **no FATAL line** — `io_wrapper.cpp` calls `MPI_Abort`
BEFORE the `std::cout` that names the file, so that message is unreachable. It is NOT an
AthenaK bug and NOT related to the physics: strace shows ROMIO doing
`statfs("bin/dhj.mhd_w_bcc.00000.bin")` on the not-yet-created file to resolve the
filesystem type, getting ENOENT on tmpfs, and failing the open despite `MPI_MODE_CREATE`.
On GPFS the same command works. **Never diagnose an AthenaK run staged on /tmp.**

## Verified working

`eos = general, general_eos = table` + `rt_ck = true` runs a clean cycle on
`inputs/mhd/deep_hot_jupiter_rt_eos.athinput` at 64x64x128, OMP_NUM_THREADS=16 on the
orion login node: both device self-tests pass (isothermal |F|/sigmaT^4 <= 9.0e-17,
transparent slab = 1), 88 column solves per cell, no source-limiter clip on cycle 1.
So the blow-up in [[dhj-general-eos-ck-blowup]] is not a start-up or table-loading failure.

Minor doc error: `PROVENANCE.md` describes `wavelengths_GCM_11.txt` as 12 descending
edges; upstream it is a count line `12` followed by ASCENDING edges. The code never reads
that file (it takes the edges from the k-table), so nothing depends on it.
