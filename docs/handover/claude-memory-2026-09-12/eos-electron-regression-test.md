---
name: eos-electron-regression-test
description: "DONE 2026-08-16 (commit b93ed468): tst/scripts/mhd/mhd_eos_electrons.py — how it pins the EOS electron fraction from both the banner and the fluid, and how to dry-run it without touching athenak/build"
metadata: 
  node_type: memory
  type: project
  originSessionId: 007bb828-7120-45a1-9fd1-74b062edff22
  modified: 2026-08-16T07:26:19.337Z
---

Written 2026-08-15/16, committed as `b93ed468`. Closes the one item CLAUDE.md blocked a
merge on: `ohmic_resistivity = eos`, `eos_metal_ionization`, `eos_metal_mh` and
`eos_metal_condensation` now have a regression test. Writing it turned up
[[eos-xe-resistivity-capped]], which had to be fixed first.

Files: `tst/scripts/mhd/mhd_eos_electrons.py`, `inputs/tests/mhd_eos_electrons.athinput`,
plus `run_output()` added to `tst/scripts/utils/athena.py` (`athena.run` pipes stdout into
a LogPipe, so the startup banner is otherwise unreadable from a test).

## The one design idea worth remembering

**Compare fluid-derived numbers against banner-derived numbers from the SAME run**, so no
electron fraction is ever hard coded. The input puts the table's grid centre exactly on the
wave's own (rho, T) — `eos_logt_min/max` chosen so `0.5*(min+max) = log10(1985.14)` — so
the banner's x_e IS the wave's x_e, and the test predicts the damping from it. Do not
retune `dens`/`pgas` without moving the temperature range with them: d ln x_e/d ln T ~ 13
here, so a 1% temperature drift is a 13% eta error.

Measured margins: absolute (eos vs constant at the predicted eta) 0.1% against a 10%
tolerance; [M/H] ratio 1.8409 vs 1.8378 against 3%. Six zero-cycle banner runs plus four
1-second fluid runs — the whole test is a few seconds.

## Dry-running it without letting run_tests.py near athenak/build

`run_tests.py` deletes and recreates `build/`, and its default (non-MPI) configure FAILS on
this branch anyway ([[branch-preexisting-breakage]]). So the test was validated against
`/orion/u/jinma/ATHENAK/build_tests` through a fake harness tree at
`/orion/u/jinma/ATHENAK/xetest/`:

```
xetest/vis -> athenak/vis ;  xetest/inputs -> athenak/inputs
xetest/tstsim/scripts -> athenak/tst/scripts
xetest/tstsim/build/src/athena -> build_tests/src/athena   # build/src must be a REAL dir
cd xetest/tstsim && OMP_NUM_THREADS=1 python3 -c "import sys;sys.path.insert(0,'.');\
  import scripts.mhd.mhd_eos_electrons as t;t.run();print(t.analyze())"
```
`build/src` has to be a real directory with only the binary symlinked, because
`athena.run` does `os.chdir` and then walks `../../../inputs/`, which a symlinked `build`
would resolve wrongly.

## Environment gotchas

- **`module load` only works inside the SAME Bash call** — `source /etc/profile.d/modules.sh;
  module load gcc/13 openmpi/4.1; <cmd>`. Without it mpicxx falls back to gcc-7 and
  `-mtune=sapphirerapids` fails. A leading `module purge` breaks the module path.
- **Do not run AthenaK from the scratchpad on /tmp**: MPI-IO gives
  `MPI_ERR_NO_SUCH_FILE` when it opens the .hst. Run on /orion. See
  [[scratchpad-not-visible-to-compute-nodes]].
- A fatal error is printed unbuffered only under `stdbuf -o0`; otherwise MPI_ABORT eats it.
- The command line can only override parameters the input file already DECLARES
  (`parameter_input.cpp:390` is fatal), hence the defaults spelled out in `<mhd>`.
- `.hst` columns: `[2]=dt`, `[11]=1-ME [12]=2-ME [13]=3-ME`. For this Alfven wave the
  perturbation is in **3-ME** (2-ME stays exactly zero).
