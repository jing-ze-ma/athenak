> **SUPERSEDED IN PART, 2026-08-25: see [[dhj-isobar-vs-shell]].** The floor conclusions
> here were measured on a constant-radius shell. On the true isobar no floor touches the
> 1e-6 bar level. The 'do NOT extend x1max' finding was nightside-dominated; the dayside
> isobar is OUTSIDE the domain and needs x1max ~ 1.36e10.

---
name: dhj-floors-for-1e-6-bar
description: "MEASURED 2026-08-25: to include the 1e-6 bar level in the dhj run, lower pfloor 1.0 -> 1e-3. x1max must NOT change and dfloor must NOT be raised. tfloor_kelvin=50 crashes the run."
metadata:
  node_type: memory
  type: project
---

Task: hydro-only, tabulated EOS + correlated-k dhj. Goal = have p = 1e-6 bar inside the
domain, unfloored, without shrinking the timestep. Work dir
`/viper/ptmp/jinma/claude_eos_gpu/fl`, binary `athena_flsw` (HEAD 3429f59f), apudev,
production grid 64 x 64 x 128, 32 blocks, ~12 min per arm (~1 rotation = 3.05e5 s).

## THE ANSWER

```
pfloor        = 1.0e-3     # was 1.0
dfloor        = 5.0e-13    # was 7.26e-12
tfloor_kelvin = 200.0      # UNCHANGED -- do not lower
x1max         = 13.04e9    # UNCHANGED
```

**1 barye = 1e-6 bar exactly** (units are cgs). So the old `pfloor = 1.0` pinned the top
of the atmosphere at precisely the level of interest. Measured in the evolved reference
run: the outermost 7 cells sit at exactly p = 1.000e-06 bar, T = 2159.9 K -- a flat
artificial slab. The old setup did not miss 1e-6 bar, it MANUFACTURED it.

With `pfloor = 1e-3` and everything else fixed: smooth monotonic profile, 1e-6 bar at
r = 1.2568e10 with 9 real cells above it, top cell 3.567e-9 bar / 200 K.
**Cost: 0.2 % more cycles per unit sim time. Free.** Mass drift identical (+0.061 % vs
+0.062 %). `efloor` (the pfloor counter) goes to EXACTLY ZERO -- pfloor stops binding
altogether; the 200 K `tfloor` takes over as the active floor, which is the intent.

## DO NOT extend x1max -- this was the trap

The INITIAL condition tops out at 7.46e-5 bar at x1max = 13.04e9, which says you need
x1max >= 14.2e9 to reach 1e-6 bar. **That is wrong.** The atmosphere expands and cools
under correlated-k, and in the EVOLVED state 1e-6 bar sits at r ~ 1.23-1.26e10, well
inside the existing domain. Always measure this on the evolved state, not the IC.

Extending anyway (13.40 / 14.20 / 14.60e9) appears to give more cells above 1e-6 bar
(17 / 26 / 28) but they are JUNK: in `x142lo` the pressure falls only 3.48e-8 -> 2.25e-8
over ten cells, a stagnant floor-supported halo. `dfloor` firings jump 6x (1.7e9 -> 1.0e10)
and mass drift flips sign (+0.06 % -> -0.17 %). Cost 1-7 % in cycles for nothing.

## Raising dfloor buys NO timestep -- do not re-propose

dfloor 5e-13 / 1e-12 / 2e-12 at x1max = 14.6e9: dt at cycle 3000 = 8.59 / 8.28 / 8.50.
No trend, all noise. The dt lever people expect here does not exist, because the IC top
is isothermal so c_s is nearly flat with height.

## tfloor_kelvin = 50 CRASHES the run (cost one full sweep)

First sweep set `tfloor_kelvin = 50` alongside the low pfloor/dfloor -- three factors at
once. All five low-floor arms died with dt exploding 10 s -> 1e9-1e19 s within ~1000-1500
cycles. Re-run at `tfloor_kelvin = 200` with everything else identical: all six clean to
3000 cycles. The table floor is 31.6 K (`eos_logt_min = 1.5`) and correlated-k cools the
tenuous top hard, so 50 K sits too near the table edge. **Keep 200 K.**

## Caveats, stated honestly

Hydro only, ~1 rotation, one replicate per arm. The arm-to-arm spread in cycles/1e5 s is
~7 %, so only the pfloor result (0.2 %, and mechanistically confirmed by efloor -> 0) is
tight; the x1max penalties sit inside the scatter and the case against extending rests on
the PROFILE being floor-dominated, not on the timing. `pfloor = 1e-3` leaves only 3.6x
margin over the evolved top pressure (3.567e-3 barye); `1e-4` would give 36x at no cost
but was not measured. For the MHD production run `max_eta = 1e13` still stands, see
[[xe-resistivity-long-runs]]; nothing here was tested with a magnetic field.

## Long confirmation run QUEUED 2026-08-25: job 11009786

`/viper/u2/jinma/ATHENAK/bench/ck_hydro_long` -- pure hydro, tabulated EOS, correlated-k,
the floors above, `x1max = 13.04e9`, `tlim = 8.64e8`, 24 h `apu1` slot, binary = HEAD
3429f59f snapshot. Pending on priority when queued (apu1 fully allocated). `NOTES.md` in
that directory carries the full rationale; `resubmit_viper.sh` continues after each
timeout from the newest `rst/*.rst`.

**Tripwire on pick-up:** `dhj.log` column `eos_efloor` must stay 0. It is the pfloor
counter, and 0 means the pressure floor is still non-binding. If it starts firing the top
has cooled into the floor -- drop `pfloor` to 1e-4 (36x margin instead of 3.6x; strictly
less binding, but not measured). Second thing to watch: mass drift in `dhj.hydro.hst`,
+0.061 % over the first rotation from `dfloor`, which cannot continue indefinitely.
