---
name: red-giant-flux-deficit-is-spinup
description: "the 0.34 L emergent flux is NOT a transport bug: grey_prod reached t=3e6 with convection spun up exactly in the layers whose turnover time is shorter than the elapsed time. Full relaxation needs t ~ 1e8 s (~40 h on 6 nodes); chain launched as prod_topre"
metadata:
  type: project
---

Measured 2026-09-08 late from `grey_prod` (job 194160, the wall run on the OLD binary
`build_rg_diag`), which **ran the full tlim = 3e6 s** at dt 30.6 flat and then hit the
known harmless Kokkos-finalize abort. Scripts in the scratchpad pattern of `topprof.py`:
`fluxprof.py` (F_enth/F_kin/F_rad vs r), `pulse.py` (<v_r>, rms, Mdot), `vmlt.py`
(v_rms against the column's own `v_mlt`). Radii come from `rtcol_t0.txt` column 1 (the
face list) -- the grid is stretched, so never assume uniform dr.

## The finding

`v_rms / v_mlt` at t = 3e6, by region, against `t_turn = H_p/v_mlt` from `column.txt`:

| region | r/R | t_turn | v_rms/v_mlt at t=3e6 | verdict |
|---|---|---|---|---|
| outer envelope | 0.93 - 1.07 | 1e5 - 1.6e6 s | **0.35 - 0.70** | spun up, converged |
| mid envelope | 0.3 - 0.7 | 7e6 - 1e7 s | 0.04 - 0.14 | never started |
| deep, on the wall | 0.11 - 0.24 | 5e6 - 8e6 s | **4 - 10, growing** | breaking out |

The elapsed time is 3e6 s. **Convection has developed in exactly the layers whose
turnover time is shorter than that, and nowhere else.** The KE growth in the hst
(KE1 2.87e42 -> 4.59e44, e-folding ~3.8e5 s, still growing at tlim) is that spin-up, not
a numerical instability: dt never moved off 30.6.

In the outer envelope `F_enth/F_req` is a flat **0.33** over i = 210-280, and the RT
emergent flux measured earlier was 0.34 L. Those are the same number: one third of L is
being carried, self-consistently, by the third of the envelope that has begun convecting.
**So the 0.34 L is a spin-up transient, not the transport failure it was read as** in
[[red-giant-open-outer-seals-the-star]]. Nothing upstream is broken.

## Two traps this cost time on

- The hst `tot-E` has only **6 significant digits**, and L*3e6 s is 6e-8 of E. Energy
  conservation CANNOT be checked from the hst on this problem; the apparent "653 L"
  heating rate is pure print quantization. Use a flux profile instead.
- `F_enth` deep down is a residual of a gross `rho v h` that is ~1e6 times larger, so it
  reads as 1800x L at i = 10-50 and means nothing there. `v_rms/v_mlt` is the diagnostic
  that actually discriminates; the flux ratio is only trustworthy in the outer envelope.

Energy input was verified correct: `rad_flux_inner = -1` sets L/(4 pi r_in^2) = 1.65e13
erg/cm2/s and `conduction.cpp:342` adds it to the inner face, matching `F_req` at i = 0.

## What it costs to finish

`t_turn` peaks at ~1e7 s around 0.4 R, so full relaxation needs **t ~ 1e8 s**, 33x what
ran. grey_prod did 3e6 s in 4423 s on 6 nodes / 96 ranks => **1474 s wall per 1e6 s of
star**, so 1e8 s is ~40 h = four chained 12 h jobs. Cutting the domain does not help:
anything below ~0.9 R has t_turn > 1e6 s.

## The chain (RUNNING at handover)

`/orion/ptmp/jinma/Athenak/red_giant/prod_topre`, restarting from
`grey_prod/rst/rg.00015.rst` (t = 3e6) on **`build_rg_prod`** (the OpenMP-ON build made
this session -- see [[orion-build-openmp-and-module-traps]]), tlim 1e8, bin dt 1e6,
rst dt 1e7, 6 nodes / 96 ranks. `problem/rt_top_re` is ON by its pgen default; do NOT
pass it on the command line, because a command-line override of a parameter that is not
in the restart's embedded input is **fatal**, not ignored (this killed job 194253).

Jobs **194259** (running) + **194262/3/4/5** chained with `--dependency=afterany`.
`sub.sh` re-restarts from the newest `rst/*.rst` it finds, so resubmitting it continues
the chain with no editing. Verified after the restart: dt flat at 30.595 (unchanged from
the wall run, so the radiative-equilibrium top did not perturb it), 0.0461 s/cycle,
mass exact, KE still climbing smoothly. Five jobs x 12 h ~ 1.5e8 s of star, enough.

**What to check first next time:** does `v_rms/v_mlt` in the MID envelope (i = 60-140)
climb off 0.04-0.14, and does the outer `F_enth/F_req` climb off 0.33 toward 1?
Those two are the whole question. Watch also the deep runaway at i = 20-40 (4-10x v_mlt
on top of the wall) -- it should saturate as those layers spin up, and if it instead
keeps growing the inner boundary needs revisiting.

Related: [[red-giant-session-2026-09-08]], [[red-giant-open-outer-seals-the-star]],
[[red-giant-grey-two-stream]], [[orion-build-openmp-and-module-traps]].

## Checkpoint 2026-09-08 23:20, job 194259 (t = 5.49e6, 2.5e6 s in, 3800 s wall)

Answer to "what to check first": **yes, the mid envelope is filling in, from the inside
out.** `tot/req = (F_enth + F_kin + F_rad)/F_req` (the kinetic flux IS in the sum -- deep
it is a large residual comparable to F_enth, in the outer envelope only -3 %, so it does
not move the 0.33):

| i | 100 | 110 | 120 | 130 | 140 | 150 |
|---|---|---|---|---|---|---|
| t=4e6 | 0.14 | 0.057 | 0.022 | 0.012 | -0.005 | 1.51 |
| t=5e6 | **2.03** | **0.73** | 0.30 | 0.10 | **0.64** | 1.04 |

The dead zone between the deep front and the surface layer narrowed from i = 90-145 to
i = 115-140 in 1e6 s, i.e. the front advances ~20-30 cells per 1e6 s. `v_rms/v_mlt` says
the same: i=60 went 2.95 -> 7.95 and i=80 0.35 -> 2.75 in that interval, while the deep
runaway at i=20-40 **saturated** (12.7 -> 12.6 at i=20), which is what it was supposed to
do. Outer envelope (i = 200-280) still flat at 0.32-0.37, unchanged -- it will only move
once the front reaches it. dt still flat at 30.36, mass exact, nothing to intervene on.

`prod_topre` has no `rtcol_t0.txt`; the scripts want one -- symlink `rtcol_late.txt` to it.
