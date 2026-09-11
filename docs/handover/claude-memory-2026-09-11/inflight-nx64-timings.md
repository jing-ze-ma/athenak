---
name: inflight-nx64-timings
description: "In-flight task as of 2026-08-16 — mesh nx2 raised to 64, GPU timings for docs sections 5/6 still to be measured"
metadata: 
  node_type: memory
  type: project
  originSessionId: 727372a4-3557-4a71-b294-5104e4fbe37f
  modified: 2026-08-19T10:50:59.814Z
---

**Task, paused mid-way on 2026-08-16:** set `<mesh>/nx2 = 64` in `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` (DONE, uncommitted) and re-measure the timings in sections 5 and 6 of `docs/general_eos_gpu_porting.md`, **on GPU**, at the new resolution. The user pointed at `run/2500_test/ohm_128/submit_viper.sh` as the submission pattern.

**Already committed and pushed** to `fork/general-eos`: `d52f3610` (host-pointer-on-device fix, 5 sites) and `0accdbdb` (meshblock table corrections). Working tree holds only the uncommitted `nx2 = 56 -> 64` edit; commit it together with the new timings.

**As of 2026-08-19:** branch fast-forwarded onto `fork/general-eos` tip `d1c289dd` (4 new commits: ideal-gas resistive input+docs, event-log MPI deadlock fix, warm-start temperature inversions). `build-gpu-bench/` rebuilt clean on top of that. Local commit `4cc3690c` (the `Resistivity::s` fix) is **committed but deliberately NOT pushed** — the user declined the push, so ask before pushing rather than assuming.

**The GPU works.** First-ever GPU run succeeded: job 10931501 on `apudev`, MI300A `gfx942:sramecc+:xnack+`, 50 cycles, exit 0, 1.43e7 zone-cycles/s. `dt_diff` matched the CPU run bit-for-bit (1.935608e+04). Section 7's "nothing has executed on an accelerator" is now out of date.

**FIXED 2026-08-19 in `4cc3690c`** (local, unpushed): `Resistivity::s` was `int s;` with no initialiser, now `int s = 0;`. Correction to the earlier note — it is assigned in `TotStage`, not `RKGCoeff`; `resistivity_rkg_tasks.cpp:83` sits inside `TotStage`, which runs in `before_rkg_timeintegrator` at `driver.cpp:410`, i.e. *before* the stage loop at `driver.cpp:412` reads `s` as its bound. So only the cycle-0 diagnostic print at `driver.cpp:554` ever read garbage (CPU `stage=0`, GPU `stage=-4096`). Not a GPU bug.

**Still to measure** (all on GPU, 64-theta, meshblock 32x32 = 8 blocks): section 5's STS table (5 tau rows x STS on/off, swept via `ohmic_resistivity = constant` with varying `eta_ohm_const`), section 5's accuracy A/B at `eta = 1e14`, and section 6's general-EOS-vs-ideal+perna cost at `bbot` = 3 G and 10 G. The existing numbers are CPU-measured at 56-theta and must be replaced or relabelled, not left implicitly current.

**How to run.** Build with [[viper-hip-build-recipe]] into `build-gpu-bench/` (gitignored by `build*/`; do NOT clobber the user's own `build/`, which is their working GPU build). **`/tmp` is node-local on Viper** — anything a SLURM job must read has to live on `/viper/u2`. Scratch run dirs go in `/viper/u2/jinma/ATHENAK/bench/`, never in [[never-write-in-run-dir]]. `apudev` is fine per the user but caps at 15 min with 2 nodes; longer rows need `apu`. Job scripts need `export HSA_XNACK=1` and `mkdir -p bin` beside the input or every rank aborts with `MPI_ERR_NO_SUCH_FILE`. Set `--mail-type=NONE` on generated jobs; the user's own script mails on ALL.

Push with the repo-scoped deploy key: `git push fork general-eos`. See [[use-fork-not-origin]].
