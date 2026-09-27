---
name: red-giant-session-2026-09-08
description: "START HERE for the red giant: what ran on 2026-09-08 (orion), the ten commits, the three jobs left running, and what to do first next time"
metadata:
  type: project
---

**Three jobs were left RUNNING at handover (12 h wall, submitted ~13:00, tlim 3e6 s):**
- `193861` taufix100 -- blend tau 10/100, `/orion/ptmp/jinma/Athenak/red_giant/taufix100`
- `193862` taufix300 -- blend tau 30/300, same but deeper. Both past t = 2.1e5, dt 30.66 s.
- `193869` combo -- **the one to look at first**: taufix100 PLUS the MLT velocity seed at
  Mach 0.05 (`vpert = 1`, `vpert_mlt = true`, `vpert_mach_max = 0.05`),
  `/orion/ptmp/jinma/Athenak/red_giant/combo`. Started clean at dt 30.66 s, no NaN.
  It is the first configuration that can BOTH survive past 1.9e5 and start with a
  realistic convective amplitude, so it is the one that can answer whether convection
  develops and whether the emergent flux relaxes toward L.

`193857` mltseed5 was CANCELLED at t = 1.83e5 with dt down to 0.039 s -- it had the old
1/10 blend. What it established before dying: **Mach 0.05 is the stable ceiling for the
MLT seed** (Mach 0.3 destroyed the run at t = 5.8e4 with a hydro blow-up spanning the
whole radial range), and its KE still DECAYED (KE1 2.87e42 -> 1.66e42, KE2 5.85e42 ->
2.85e40), i.e. a bigger seed alone did not start convection.

**READ FIRST:** [[red-giant-dt-collapse-solved]] -- the control set and the fix.

## The production configuration that works
cubed sphere, `-DPROBLEM=red_giant`, build `build_rg_diag/`. 6 x 32 x 32 x nx1=320,
meshblock 320 x 8 x 8 = 96 blocks. Wall at r_in = 9.6e10 = 0.03 R (BELOW the convection
zone: the RCB is at 0.04 R, measured), `rad_flux_inner = -1`, MLT column at
`mlt_alpha_ic = 3`, monotone stretch (`--dumin 0.2`), EOS table to logT 7.5 / logrho 2,
`ic_tau_rad = 10`, blend `rad_tau_lo/hi = 10/100`, sponge on (`sponge_zbot = 0.98`).

## SLURM on orion (this cost real time)
`--cpus-per-task` counts **CORES (112/node), not the 224 hyperthreads** SLURM reports.
16 ranks x 7 = 112 fills a node; asking 24x9 or 16x14 is rejected as "Requested node
configuration is not available". Measured at 96 ranks: 4 nodes (96/112 cores)
0.0747 s/cycle, **6 nodes (full) 0.0547**, 8 nodes 0.0625. 96 ranks is the ceiling --
the radial direction CANNOT be split because the two-stream sweeps whole columns.

## Ten commits, in order
979b7f61 conduction dt names the cell AND its state; open-BC energy budget + `open_conserve`
91514722 `problem/face_budget` -- what the radial faces advect; it found the drain
ce9130b6 open BC: zero the NET mass flux through the inner FACE (the cell-centre pass missed it)
7230fd53 the same pass ported to solar_convection (it gains 4.7 % of its mass in 1.2e4 s)
c0634e7e `mlt_alpha_ic` -- build the column on the MIXING-LENGTH gradient, not grad_ad
6ef9492d **the first cycle ran with NO conduction timestep** (blend weights not built yet)
7c2a2268 `wall_noflux` (the wall was NOT impermeable, -126 L) + a top sponge
19a5f645 `fit_radial_stretch.py --dumin` -- a MONOTONE fit; plain LSQ folds over 22 scale heights
d9992676 MLT-amplitude velocity seed (`vpert_mlt`, `vpert_mach_max`)
25692d16 `ic_tau_rad` -- the initial column no longer depends on the RT handover

## What to do first next time -- SPENT, see [[red-giant-flux-deficit-is-spinup]]

The four items below were all worked on 2026-09-08 night. 193861/2/869 are finished or
cancelled; `grey_prod` superseded them and reached tlim 3e6. Item 3 (promote the working
input to the tracked inputs) is the ONE still open. The current production input is
`/orion/ptmp/jinma/Athenak/red_giant/prod_topre/rg.athinput`.

## The original list
1. Check 193861/193862: did they reach t ~ 1e6, and does `tools/grid/flux_gate.py`
   (on `rtcol_t0.txt` / `rtcol_late.txt`) show the emergent flux moving TOWARD L?
   0.877 -> 0.556 was the old failure and it is still the open physics question.
2. Does convection actually develop -- is KE growing or still decaying?
3. Nothing is committed to the tracked INPUTS yet; the working input is
   `/orion/ptmp/jinma/Athenak/red_giant/taufix100/rg.athinput`. Promote it.
4. The MLT velocity seed blows the run up above ~Mach 0.05 with a single smooth mode
   (`mltseed4` died at 5.8e4 at Mach 0.3). If it is wanted, it needs a multi-mode seed.
