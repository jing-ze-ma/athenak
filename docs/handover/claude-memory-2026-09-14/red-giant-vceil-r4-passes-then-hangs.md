---
name: red-giant-vceil-r4-passes-then-hangs
description: "VCEIL DONE 09-11 07:00: <hydro>/vceil (Newtonian, default 0; c2p off the cs, metric-norm deferred path in GnomonicEquiangleRaiseVel on the cs; counter eos_vceil) + RT NaN guard in eiN (rt_use_cons): pin10. Gates pass (off bit-identical to pin9; healthy star active max |v| 2.5e6 << 5e7, counter fires only in GHOSTS). R4_vceil (5e7) PASSED the 5.87e5 death with dt flat 30.63 to 6.13e5, then HUNG in MPI at cycle 20000-20002 (no error); I8_vceil 196221 launched to 9e5 (may hang at 20000 too); hang-diagnosis agent running"
metadata:
  type: project
---
Diff verified by me: fs = vceil/|v|, m *= fs, u.e -= (1-fs^2) KE, e untouched; cs path uses
|v|^2 = 2 ekin/d with the metric ekin and writes u0(IM*) immediately; skipped on only_testfloors.
R4 numbers (t, dt, collapses, crossings, eos_vceil/window, L_out/L): 5.87e5 30.63 0 ~190 5321 ~150;
6.0046e5 30.63 0 205 9519 4.16; 6.1266e5 30.63 0 230 10607 2.15. RT clamp count 0 (the NaN path
never triggered). The ceiling fires ~4/cycle on the healthy star but only in ghost cells ->
hygiene fix queued (restrict to active cells). Runs: R4_vceil (196207), I8_vceil (196221),
gates GV1a/GV1b/GV2a/GV2b. Binary athena_pin10 (needs the SIX rt opt-ins in the input).
Physics next (user agreed 07:30): raise the corona base density 1e3-1e4 (to ~1e-14 at the join,
T 2e4 K stays hydrostatic), raise rad_gate_rho just above it, keep vceil as guard; the shell
launch at 5.7e5 itself is still unexplained. See [[red-giant-r3-dfloor-keep-velocity-not-the-loop]].

## HANG RESOLVED 07:40: runaway_scan.hpp line budget was rank-local (owner rank alone
## incremented `lines`; at 2000 it returned before the MPI_Allreduce(MAXLOC) -> deadlock, like
## eventlog 27da6380). Fixed (++lines before the owner filter). pin11 = pin10 + that + vceil
## restricted to active cells (eos_vceil 0 on the healthy star, dumps bit-identical, GV3a/b).
## I8_vceil (no scan) never hung: 6.26e5 dt 30.63 at 07:40.

## C1_dense_corona LAUNCHED 08:25 (jobs 196255 + chain 196256): from scratch, pin11, I8 input +
## bg_rho 1e-15 -> 1e-12 (bg_rho is the JOIN threshold; rho_corona = 0.081 bg_rho -> 8.43e-14 at
## the join, 1041x), rad_gate_rho 1e-16 -> 8.1e-13, rst dt 1e5. CAVEATS: the join moved inward
## 3.652e12 -> 3.556e12 (corona replaces the outermost ~3% of the atmosphere; last stellar cell
## 1e-12, 3364 K); the gate ramp partially gates the top 2-3 stellar cells (G 0.69 at the last).
## 2.74e5 sim s per wall hour -> 9e5 in ~3.3 h. Compare vs I8: eos_vceil counts, shell at 5.7e5.

## I8_vceil REACHED 9e5 (09:40): first open-top run of the lineage to finish; dt 30.63 flat
## except 7 one-cycle dips 6.03-6.25e5 (dfloor cells rho=1e-18, v=0, T 1e16: the DENSITY floor
## keeps the ENERGY -> T = e/rho explodes; hygiene fix = scale e with rho, "dfloor_keep_temperature");
## L_rad,out/L 2.89 at 8.08e5, 1.93 at 9e5 (I5: 3.3 at 8.2e5, 2.0 at 9e5); eos_vceil ~9500/window
## at the end. C1_dense_corona at 2.88e5: 0 collapses, eos_vceil 0 so far, L_out 0.66.
## Post-run analysis agent -> _analysis_0911/i8_report.txt.

## I8 ANALYSIS (09:55, _analysis_0911/i8_report.txt): HEALTHY at 9e5. tot-E/mass constant to hst
## precision (<1e-5) incl. the shell interval; L_out/L 1.74/3.30/1.93 at 7/8.25/9e5 (I5 1.2/3.3/2.0),
## lidded band 1.1-1.3 between episodes; shell BREATHES (r(1e-12) 3.59->3.70e12 by 8e5, back to
## 3.66e12 at 9e5), not escaping. Ceiling: 1.58e6 firings total, peak 2e5/interval at 6.6e5, but NO
## dump has |v| within 5x of 5e7 (max 1.07e7) -> clips within-cycle transients in the corona above
## the front (r>3.75e12), guard only. OPEN BUGS: (a) 348 RT cell reports 6.67-7.13e5 with -nan
## I_dn/F/src in OPTICALLY THICK columns (tau_to_top 75-90), all caught by the rescue (de/e=-0.999)
## -> the NaN is not the T guard, chase it in the thick-limit two-stream algebra; (b) the 7 dt dips
## = dfloor cells keeping ENERGY (T 1e16) -> scale e with rho on the density floor.
