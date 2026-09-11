> **RETIRED 2026-08-25, superseded by [[ck-limb-run]].** Job 11009786 stopped at
> t = 1.184e7 s (137 days). Its grid could not put the 1e-6 bar isobar inside the domain
> at the terminators. The floor numbers below were measured on a constant-radius shell,
> which is not the 1e-6 bar level -- see [[dhj-isobar-vs-shell]].

---
name: ck-hydro-long-run
description: The ck_hydro_long confirmation run -- 1000-day target, ~67 h wall, 3 slots; mass drift is a transient not a leak
metadata:
  type: project
---

`/viper/u2/jinma/ATHENAK/bench/ck_hydro_long`, job 11009786 on apu1, started 2026-08-25
13:28. Pure hydro, tabulated EOS + correlated-k, binary = HEAD 3429f59f. Confirms the
floor setting from [[dhj-floors-for-1e-6-bar]] over many rotations.

**Simulation "days" are Earth days (86400 s), not rotations.** The rotation period is
P = 3.05e5 s = 3.53 days. The athinput `tlim = 8.64e8` is exactly 10000 days -- 10x more
than needed and ~28 days of wall. The user's targets: **500 days is analyzable, 1000 days
is the better treatment.**

Measured slot-1 throughput: **357 sim-s per wall-s**, stable once dt settles at ~6.2 s
after rotation ~6 (starts at 13.4). So 500 days = ~34 h total (2 slots), 1000 days =
~67 h total (3 slots), i.e. only two resubmits -- not worth automating.

`tlim` is capped at 1000 days via a **`time/tlim=8.64e7` command-line override** on the
`srun` line of `resubmit_viper.sh`, not by editing the athinput. Verified safe on restart:
`main.cpp:275` `ModifyFromCmdline` runs AFTER `main.cpp:263` `LoadFromFile`, and
`driver.cpp:88` reads tlim from `pinput`.

**Two NOTES.md watch items resolved in slot 1:**

* `eos_efloor` stayed **0** through 841k cycles -- `pfloor = 1e-3` is non-binding as
  measured, no need to drop to 1e-4.
* **Mass drift is a transient, not a leak, and not per-rotation.** It rises +0.062 %
  during rotation 1, plateaus at +0.078 % by rotation 10, then DECREASES (+0.073 % at
  rotation 16.7). The old worry "+0.061 % per rotation, not sustainable" was reading a
  one-off startup transient as a rate.

**Flow checked at 50 days (dump 00002), healthy.** Equatorial superrotating jet formed,
peak zonal-mean v_phi **+2.34 km/s** at lat -2.3 deg, saturated aloft (2.296 -> 2.337 km/s
between 25 and 50 d) while total rho*v^2 fell 15 %. No meshblock imprint, polar rows stable
(polarmost/next-row max|v_theta| = 0.97, unchanged), all fields finite. The deep atmosphere
below r ~ 1.1e10 is still at rest -- an independent argument for 1000 days over 500.

**Pressure recovered via [[eos-table-dump]].** Domain spans **66 bar to 10^-6.6 bar**; the
1e-6 bar level is at **r = 1.266e10 cm**; the jet core sits at 10^-5.0 bar, a decade
deeper. T = 1648 K at 1 bar; at 1e-6 bar it is **1677 K dayside vs 856 K nightside**.

**Caveat: dfloor is active at the observing level.** It first bites at 10^-3.5 bar and holds
**21 % of all cells at 1e-6 bar** (0 % dayside at every radius, ~2x the 21 % on the night
side). So dayside spectra at 1e-6 bar are clean; the night side and hence the day-night
contrast are density-floor-supported. The earlier sweep only tested RAISING dfloor;
lowering it is untested and costs dt (less density at fixed p -> higher c_s).

Artifact with both figures: https://claude.ai/code/artifact/6089bbdb-acc2-443e-a40f-ec2fb15d84fc
