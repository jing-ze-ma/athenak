---
name: session-state-2026-08-24
description: START HERE - the two ckrepro reproducers ran 10 h clean, no blow-up; dt plateaus at ~2 s and max_eta is irrelevant; brute force cannot reach the failure
metadata:
  type: project
---

Entry point for the next session. Branch `polar-average-perf` at `c1d513c2`.
Predecessor: [[session-state-2026-08-18]].

## Reproducer result (read at 2026-08-24 ~10:05, ~9.8 h into a 12 h slot)

| job | dir under `/orion/ptmp/jinma/Athenak/` | max_eta | cycle | time | dt |
| --- | --- | --- | --- | --- | --- |
| 190881 `dhjck_gen` | `ckrepro` | 1e14 | 308800 | 6.887e5 | 2.03 |
| 190882 `dhjck_e13` | `ckrepro_eta13` | 1e13 | 305900 | 6.484e5 | 1.41 |

**Neither blew up. No NaN, no FATAL.** Both are numerically stable and boring:
- dt falls 13.4 -> ~2.5 over the first ~8000 cycles, then **plateaus at ~2 s and stays
  there for 300 000 cycles**. Not a runaway; a new equilibrium.
- **max_eta makes no difference to dt.** The earlier "1e13 runs at 2x the dt" note was
  read off `dt_diff` only; the actual `dt` trajectories of the two runs are
  indistinguishable. After ~cycle 8000 the limiter is the hydro CFL, not the resistive
  one (`dt_diff` 1.58 vs `dt` 2.03 at 1e14).
- **The `rt_de_max` source limiter DOES fire**, in both, at cycle ~600 / t = 8.0e3, in
  112 cells -- and only once (it warns once). This answers one of the open questions;
  it fires identically with and without the blow-up, so it is not diagnostic on its own.

## Consequence: brute force is dead

Throughput is 19.5 simulated s per wall s at 16x7. `tlim = 8.64e7` would take **~51 days
of wall time**. 10 h bought t = 6.9e5, under 1 % of tlim. Chasing the viper NaN by just
running longer is not viable, so the next step must be one of:
1. **Recover the viper failure time** (cycle / `time=`) from the user's viper logs -- not
   present anywhere on orion; `run/dhj_pole_rt_viper*/` holds only the athinput. This is
   still the single most valuable missing number.
2. A cheaper reproducer: lower resolution / smaller radial range, so t advances faster.
3. Restart-chain the existing runs from `rst/dhj.00001.rst` (t = 4.32e5) in further 12 h
   slots -- costly, and blind without (1).

## The main measured result (unchanged)

**The blow-up is NOT in the initial state, and the table EOS is the SMOOTHER of the two
at t = 0** -- max neighbour |dT|/T 0.026 (general) vs 0.067 (ideal), max T in the CK
region 4315 K vs 4794 K, Gamma_1 smooth. Detail and the dead hypotheses:
[[dhj-general-eos-ck-blowup]].

## Uncommitted work in the repo (3 modified files, nothing committed)

- `inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput` -- ADDED the whole `rt_ck` / `ck_*`
  block (overrides only modify parameters that already exist, so CK could not be switched
  on for the ideal-vs-table comparison at all).
- `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` -- added an `<output3> file_type = rst`
  block (same reason), rewrote the stale `<meshblock>` comment, added the measured
  ranks x threads table.
- `data/exo_fms_ck/PROVENANCE.md` -- upstream defects and the clean checks.

CK tables installed and md5-verified in `data/exo_fms_ck/` (gitignored): [[exo-fms-ck-tables]].

## Also settled

- **16x7 is measured-optimal**; a 64^3 meshblock silently caps a run at 2 ranks.
  [[freya-job-submission]].
- Output goes to `/orion/ptmp/jinma/Athenak/`: [[test-output-location]].
- Quota cleanup freed ~7 GB on `/orion/u`: [[gpfs-quota-wall]].
- The one upstream CK typo (`cia/H2-_ff.txt`) sits in a file the code never opens.
