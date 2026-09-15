---
name: session-state-2026-08-16
description: "Where everything stood at the end of 2026-08-16: four unpushed commits on general-eos, three threads closed, one still open (the solar atmospheric runaway). START HERE."
metadata:
  type: project
---

End of 2026-08-16. **Read this first, then the linked note for whichever thread is being
picked up.** Working tree clean under `src/ inputs/ tst/ docs/ tools/`.

## ALL PUSHED — `origin/general-eos` is at `5d1c3436`

Verified with `git ls-remote`. Working tree clean under
`src/ inputs/ tst/ docs/ tools/`. Upstream untouched, no PR.

Nothing is outstanding. Ten commits on 2026-08-16, oldest first:

- `c0edbe49` fix two out-of-bounds reads of `eta_b` (one made `ohmic_resistivity = constant`
  a silent no-op) — [[eos-xe-resistivity-capped]]
- `b93ed468` regression test for the EOS electron fraction — [[eos-electron-regression-test]]
- `ae33271d` solar_convection tools: tag acoustic/lidtest outputs by run set
- `cabf4922` solar_convection: `problem/z_ph_frac` — where the photosphere sits in the box
- `c8657f34` `tfloor_kelvin` — the temperature floor in kelvin ([[general-eos-project]])
- `bc6b5774` deep_hot_jupiter_rt: drop the outer-BC Maxwell term
- `767fff43` deep_hot_jupiter_rt: put it back as a clamped effective gravity (supersedes it)
- `9c06ea72` reference dhj EOS input carries the validated settings (max_eta 1e14 + STS)
- `ee968145` `docs/general_eos_gpu_porting.md`
- `2dedcbdb` general EOS: drop redundant table interpolations (10%, bitwise-verified)
- `5d1c3436` ideal gas can use the EOS electron fraction; 7% FASTER than perna
  ([[resistivity-xe-table-for-ideal]]) -- also fixed an uninitialised `HotJupiterParam`
  whose Rgas the resistivity divides by, and a `GetOrAddString` probe that switched the
  resistivity on for runs that had none

Working tree clean under `src/ inputs/ tst/ docs/ tools/`. Upstream still untouched, no PR.

## Thread status

| thread | state | note |
|---|---|---|
| EOS electron fraction / `ohmic_resistivity = eos` | **CLOSED** | [[eos-xe-resistivity-capped]], [[eos-electron-regression-test]] |
| dhj at high B, diffusive dt, STS | **CLOSED** — BC reworked properly, STS A/B done, general-EOS cost measured (~3x) | [[dhj-highB-outer-bc]] |
| `tfloor_kelvin` | **CLOSED** | [[general-eos-project]] |
| solar_convection atmospheric runaway | **OPEN** | [[solar-convection-general-eos]] |
| general-EOS cost / Stage 4 (rho,e) table | **DECISION DEFERRED** by the user 2026-08-16 — scoped, not started, do not begin without asking | [[general-eos-stage4-rho-e-table]], [[general-eos-table-cost]] |

## The one open question

The solar atmosphere does not reach a steady state. Corrected sponge: negative. Taller box
(16.7 scale heights instead of 8.7): delays the runaway to a quasi-steady t~2000-15000 and
then it resumes; dE/dt/F+ is still positive in all four runs at t=20000. **Next candidate
is a real energy sink — radiative damping in the atmosphere** (t_rad/t_ac is 20-45 there,
i.e. undamped). Before quoting the tall runs, resolve why their atmosphere is ~1.9x less
dense than the short runs' at every height above the photosphere.

## Scratch directories on /orion (none are in [[run-directory-untouchable]])

- `soleos/` — all solar_convection runs; `tall/{nosponge,sponge}` are the new ones,
  `sponge2/` the corrected-sponge set. `soleos/tools/*` symlink to `athenak/tools/`.
- `dhjb/` (~420 MB; the .bin dumps were deleted, logs / inputs / md5 records kept) — the whole dhj
  high-B investigation: `stress/` (the bbot x dfloor sweep that
  reproduced the crash), `loc/` (dumps that localised it to the last two radial shells),
  `vac/` (the reverted Alfven-ceiling test), `nomax/` (the fix, verified), `bench/` (the
  STS cost benchmark), `b3cmp/` (3 G on-vs-off comparison), `prof/` + `perfrun/` (the cost
  decomposition and the perf profile), `tol/` (logtol and no-root-find ceiling), `xe/`
  (ideal-gas x_e table), `bw/` (bitwise md5 records, dumps deleted). Binaries `athena4`
  (before the EOS optimisation) and `athena6` (current) are kept for A/B; the rest were
  deleted.
- `xetest/` — the Alfven-wave resistivity tests AND `xetest/tstsim/`, the fake test-harness
  tree that runs `tst/` modules against `build_tests` without letting `run_tests.py` delete
  `athenak/build`. Recipe in [[eos-electron-regression-test]].
- `build_tests/` (built-in pgens), `build_dhjrt/` (PROBLEM=deep_hot_jupiter_rt),
  `athenak/build/` (PROBLEM=solar_convection). All three are current as of this date.

## Cluster habits confirmed today

- `module load` only works inside the SAME Bash call: `source /etc/profile.d/modules.sh;
  module load gcc/13 openmpi/4.1 cmake/4.0; <cmd>`. A leading `module purge` breaks it.
- **AthenaK does not create its output directory**: `mkdir -p bin` in the run dir first or
  every rank dies with `MPI_ERR_NO_SUCH_FILE`.
- p.shared with `--ntasks-per-node=16 --cpus-per-task=7` uses all 112 physical cores of a
  node and divides 128 meshblocks evenly. Starts within a minute or two.


## Running this on Viper (2026-08-16)

**MEMORY DOES NOT TRANSFER TO VIPER.** `/u/jinma` is the ORION filesystem (`df` -> `orion_u`
on `/orion/u`), which Viper does not mount; and the memory dir is keyed by the project's
absolute path (`-orion-u-jinma-ATHENAK-athenak`), so a repo at any other path gets a
different, empty one. The Claude on Viper starts cold.

**Git is the handover channel, and `run/` is untracked** — so the validated viper input in
`run/dhj_pole_rt_viper/` does NOT travel. That is why `9c06ea72` put the settings into the
tracked `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` and `ee968145` added
`docs/general_eos_gpu_porting.md`. Between them they carry the build flags
(Kokkos 4.6.2 knows `AMD_GFX942`), the meshblock advice, every measured number, and the
GPU-portability defect. Point a fresh agent at those two files.

**The defect worth remembering:** `resistivity.cpp` captures `this` in ~a dozen device
lambdas and `CurrentDensity()` in `current_density.hpp` takes `MeshBlockPack*` and reads
`pmesh`/`pcoord` through it inside a `KOKKOS_INLINE_FUNCTION`. Legal on MI300A (APU,
hardware-coherent memory) — which is the ONLY reason it runs there — and illegal on any
discrete GPU. Afflicts `perna` as much as `eos`. Fix is contained: `current_density.hpp`
has exactly one includer.
