---
name: drop-unhelpful-defaults-0926
description: User 2026-09-26 -- default-on switches measured as NOT helping (on H200) should be dropped from default-on; candidate so far only M1 implicit_one_pass (pending multi-config measurement); may need a backend-dependent default (it helped on MI300A)
metadata:
  type: feedback
---

User, 2026-09-26: "for the switches that we measured are not helping, drop them from default on".

**Status (09-26 morning)**
- dhj ck levers all help on H200. cvsec is cost-neutral but needed (off gives 51/300 non-converged).
- mg_gf is dropped.
- **The only candidate is `<rad_m1>/implicit_one_pass`:** -3.5 % when off on the He box, 1 H200. The one_pass agent is measuring more configurations (work dir /resnick/scratch/jingze/onepass/).

**How to apply**
- Flip a default off only on multi-configuration evidence, with the usual gates:
  - old value explicit = bitwise;
  - key absent = new value;
  - other problems bitwise.
- If a switch helps on MI300A but hurts on H200, propose a backend-dependent default (`KOKKOS_ENABLE_CUDA`) instead of removing it for viper.
- Record each drop in the handover.

**one_pass measured 09-26 (good node hpc-sm-02-10, 150 cycles, fresh IC).**

What it does: when the predictor is on, the first Picard pass runs at the full linear tolerance, and the solve may be accepted after that one pass if `res0*qe/(1-qe) < tol`, with safety 30 for hesdirk2. So it costs extra Krylov work on every solve and saves second passes only when solves are accepted.

| config | GPUs | change | accuracy (x round-off) |
|---|---|---|---|
| He box, cfl 0.3 | 1 | **-5.5 %** | |
| He box, cfl 0.3 | 2 | **-10.6 %** | |
| reduced box | 1 | -4.4 % | |
| He box, cfl 0.9 | 1 | +4.6 % (0 solves accepted) | |
| sph_atm | 1, 2 | ~0 | ~1x |
| gateC wedge | 1 | ~-2 %, noise (0 accepted) | |
| box / gateC (default vs off) | | | **5-22x the round-off spread** (~1e-10..1e-9 relative, below implicit_tol 1e-8) |

- On speed, one_pass helps at the reference cfl 0.3, so do NOT drop it for speed.
- On accuracy, it breaches the user's round-off bar. User decision pending: keep it (tolerance-level) or set 0 (costs 5-11 % on the box).
- Possible code change: auto-disable one_pass when acceptances stay at 0 over a window.

**MERGED 09-26: implicit_one_pass_auto (default true), dff69bc7.** Where no solve is accepted after one pass, one_pass switches off per solve kind; it re-probes every 64 solves. It is bitwise wherever it never triggers. The H200 timing job 3526294 is pending.

**H200 timing of the auto-disable** (hpc-sm-02-09, 1 H200, 2 reps):

| config | old default | one_pass 0 | auto |
|---|---|---|---|
| He box cfl 0.9 | ~72.7 (one rep 83.9 outlier) | ~67 | ~69 (+2 % vs off: first window + probes) |
| He box cfl 0.3 | 50.5 | 53.9 | 50.2 (bitwise vs default) |
| gateC | 12.0 | 12.2 | 12.2 |

- gateC: the old default was NOT slower here. The tight first pass cut Picard passes 3.00 -> 2.67 even with 0 solves accepted, so auto gains nothing on gateC.
- The one_pass restart header grows 9 -> 18 Reals when auto is on.
