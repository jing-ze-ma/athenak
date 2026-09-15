---
name: ck-limb-run
description: THE LIVE CAMPAIGN -- ck_limb, sized so both terminators reach 1e-6 bar with 3 scale heights of buffer
metadata:
  type: project
---

`/viper/u2/jinma/ATHENAK/bench/ck_limb`, **job 11022482 on apu1**, launched 2026-08-25.
Binary = the SAME `athena` snapshot as the retired run (HEAD 3429f59f), so only the input
differs. **Replaces [[ck-hydro-long-run]]**, which was stopped at t = 1.184e7 s (137 days).

**The science goal, as stated by the user:** the leading and trailing terminators resolved
to the 1e-6 bar level, with the observable surface unaffected by the outer boundary. The
old grid met neither -- see [[dhj-isobar-vs-shell]] for the measurement.

**Changes from ck_hydro_long, with the reason for each:**

| parameter | old | new | why |
|---|---|---|---|
| `x1max` | 1.304e10 | **1.433e10** | 3 scale heights above the 99th-pct terminator isobar |
| `nx1` (mesh AND meshblock) | 64 | **88** | holds dr at 5.557e7; extending at fixed nx1 coarsens the top 32 % |
| `dfloor` | 5e-13 | **5e-14** | measured free, 0.0 % dt cost even at 100x (jobs 11019549-51) |
| `output3/dt` | 2.16e6 | **3.05e5** | one dump per rotation, so the limb can be time-averaged |
| `tlim` | 8.64e8 | **8.64e7** | 1000 days; days are EARTH days (86400 s) |

| `f_stretch_theta` | 2.0 | **3.0** | **MEASURED 1.407x** on apudev (11022124 vs 11022332) |

Deliberately NOT changed: `pfloor = 1e-3`, `tfloor_kelvin = 200`.

**The f_stretch_theta win is real and measured, not estimated.** dt is set by the AZIMUTHAL
CFL at the pole (`dx3 = r sin(theta) dphi`), so widening the polar-most cell relaxes the
binding constraint directly: mean dt 8.934 -> 12.558 s with cycles/s identical at 44.3 (same
cells, same work; the whole gain is dt). It also refines the equator 1.55 -> 0.84 deg, where
the terminators are, and the x2 CFL has a 10x margin so that is free. Stopped at 3.0 rather
than 4.0 because 4.0 makes the equatorial cell 7:1 anisotropic against dphi = 2.81 deg and
the pole:equator ratio 26x. **x3 is still the limiter (x1 = 59.1 s, x2 = 219, x3 = 37.1), so
there is still a 1.6x margin if a future run needs it** -- but re-check the limiter first.

**Validated on apudev before launch (job 11022124, 0.2 rotation from the IC):**
throughput 44.3 cycles/s vs the old grid's 57.9 = 0.765x, matching the 1.375x cell count,
so **dt per cycle is unaffected**; with the stretching folded in, expect roughly **1.08x the
old run's throughput** and **~3 slots to 1000 days** (confirm once dt settles). Terminator columns with the isobar above the top: **0.0 %**, vs 43 % on the old
grid; 17 cells of headroom at BOTH limbs (~3 scale heights, the design target).
The IC is the worst case: the atmosphere cools and contracts under correlated-k, which
drops the top-cell pressure further (the old run went 7.5e-5 -> 3.8e-6 bar dayside).

**Watch:** cells above the 1e-6 bar isobar at lon = +-90 (must never reach 0); whether the
new top stratifies or forms a stagnant halo (**if it does, check `eos_dfloor` first** --
the whole reason the earlier x1max test failed was that the floor pinned that halo);
`eos_efloor` is NON-ZERO here and that is EXPECTED -- lowering dfloor shifts floor
activity from the density floor to the pressure floor (A/B counters: efloor 0 -> 15.4M ->
36.7M as dfloor went 5e-13 -> 5e-14 -> 5e-15). The test is WHERE it fires: at 5e-14 every
p-floored cell sat at log10 p = -9.4..-9.0 bar and 0.00 % of cells at the 1e-6 bar level
were p-floored. That 0.00 % at the terminators is the number to watch. The substellar column is the tightest at 3 cells of headroom in
the IC, but the goal is the limb, not the substellar point.

## Verified after 1 rotation (2026-08-25)

Goal met: **19 cells above the 1e-6 bar isobar at BOTH terminators**, 0.00 % of terminator
columns above the top, and **0.00 % of cells ON the isobar at either floor**. eos_fail = 0.

**Partly falsified my own halo hypothesis, keep this.** I claimed holding dr fixed plus a
lower dfloor would let the extended top stratify. Measured H_observed / H_hydrostatic at the
top: **dayside 1.3** (fine -- genuinely hot, 1905 K, mu = 1.51 from H2 dissociation),
**nightside 11.2** with 85 % of the top ten cells at dfloor. The nightside halo is NOT a
tunable floor artifact: its true scale height there is **5.12e7 cm = 0.92 of a cell** at
dr = 5.56e7, so the grid cannot represent it at any floor value. Lowering dfloor 10x only
moved the ten-cell falloff from 1.5x to 2.4x.

**Structural tension to remember:** the dayside needs x1max = 1.433e10 to reach 1e-6 bar;
at that radius the nightside is 589 K with H < 1 cell and can only be floor-supported. One
spherical grid with a single x1max cannot satisfy both. Acceptable only because the halo is
three decades above the observable and never touches the isobar -- verify that, not the
floor counters.

Mass drift +0.117 % over 1.1 rot, ~1.9x the retired run (taller domain, 10.29 % of cells at
dfloor). The old one plateaued by rotation 10; if this does not, revisit.
