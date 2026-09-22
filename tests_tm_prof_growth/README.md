# tm_prof_growth — profile of the prod4 production form, and the split-lag growth gate

Binary for all three arms: `bench/cksd_0922/src_snapshot/build_gpu/src/athena`
(content of `rt-integration` `5f348164`), md5 `a6dd8e0bab4385972b453b03975dadd4`,
`strings -a | grep -c ck_sweep_form` = 8. 2 ranks / 2 GPUs on apudev, 24 MeshBlocks,
12 per GPU, as production.

**Binary caveat, applies to every absolute number below.** This snapshot predates both
the commit that flipped `ck_sweep_form` to `tm` by default and the FOP-tag commit, so
(a) `ck_sweep_form = 1` had to be written explicitly into the inputs — harmless, the run
logs confirm it took (`correlated-k PROBE-FREE sweep ON (problem/ck_sweep_form = 1, tm:
transfer-matrix/adding in (u,d))`) — and (b) the absolute cycles/s are about **9 % below
what HEAD would give**. Percentages, shares and the arm-to-arm ratios are unaffected;
only the headline rates should be scaled up by ~1.09 before quoting them as HEAD numbers.

---

## (7) Profile of the production (prod4) form — `prof_tm`, job 11941316

Arm `prof_tm`: restart `bench/cs_mhd_prod3/rst/dhj.00567.rst` (ncycle 4429831, rot 283),
`time/nlim = 4430131` = 300 cycles = 600 RK stages, run under
`rocprofv3 --kernel-trace --stats --output-format csv` exactly as
`bench/prof_0922/prof/submit.sh`. Input = `bench/prof_0922/plain/deep_hot_jupiter.athinput`
with the trailing `<problem>` block turned into the prod4 form: `ck_spherical = true`
(unchanged), `ck_beam_sph` false → **true**, `ck_sweep_form = 1` added.
`rt_cell_report` / `rt_report_every` stay off. All three switches confirmed in `prof.log`:
pseudo-spherical stellar beam ON, correlated-k SPHERICAL form ON, PROBE-FREE tm sweep ON,
RT split path ON with 22 chain blocks of 4. `MeshBlock-cycles = 7200` = 24 × 300, so the
arm ran exactly the intended 300 cycles.

### Top-15 kernels, rank 0

Total GPU kernel time over all **96** distinct kernels, rank 0: **18855.65 ms**
(70091 launches).

| # | % GPU | ms | calls | kernel | file |
|---|---|---|---|---|---|
| 1 | 70.34 | 13262.9 | 600 | `two_stream_rt::picket_fence_two_stream_RT_pass` (dominant specialization — the tm column sweep) | src/utils/two_stream_rt.hpp |
| 2 | 3.16 | 596.7 | 601 | `Conduction::BuildRadWeights` | src/diffusion/conduction.cpp |
| 3 | 2.53 | 477.4 | 600 | `MeshBoundaryValuesFC::PackAndSendFC` | src/bvals/bvals_fc.cpp |
| 4 | 2.30 | 433.6 | 600 | `two_stream_rt` `rt_apply` (`par_reduce_clip4`) | src/utils/two_stream_rt.hpp |
| 5 | 2.01 | 379.9 | 600 | `Conduction::ImplicitRadialUpdate` | src/utils/two_stream_column_ck.hpp |
| 6 | 1.80 | 338.5 | 601 | `GeneralMHD::ConsToPrim` | src/eos/general_mhd_frozen.cpp |
| 7 | 1.75 | 329.4 | 600 | `mhd::MHD::CalculateFluxes<hlld>` (variant 1) | src/mhd/mhd_fluxes.cpp |
| 8 | 1.64 | 309.5 | 601 | `Coordinates::GnomonicEquiangleRaiseVelMHD` | src/coordinates/gnomonic_raisevel.hpp |
| 9 | 1.60 | 301.7 | 600 | `mhd::MHD::CalculateFluxes<hlld>` (variant 2) | src/mhd/mhd_fluxes.cpp |
| 10 | 1.45 | 273.6 | 600 | `SourceFunc` (pgen source terms) | src/pgen/deep_hot_jupiter_rt.cpp |
| 11 | 1.42 | 268.2 | 600 | `MeshBoundaryValuesCC::PackAndSendCC` | src/bvals/bvals_cc.cpp |
| 12 | 1.16 | 218.7 | 600 | `mhd::MHD::CalculateFluxes<hlld>` (variant 3) | src/mhd/mhd_fluxes.cpp |
| 13 | 0.87 | 164.0 | 600 | `picket_fence_two_stream_RT_pass` (2nd specialization) | src/utils/two_stream_rt.hpp |
| 14 | 0.60 | 112.9 | 600 | `Conduction::ImplicitRadialUpdate` (2nd specialization) | src/utils/two_stream_column_ck.hpp |
| 15 | 0.58 | 108.9 | 600 | `picket_fence_two_stream_RT_pass` (3rd specialization) | src/utils/two_stream_rt.hpp |

As in `prof_0922`, the mangled names distinguish the enclosing function but not the
per-call-site `par_for` string label, so `rt_pre`, `rt_chain_ck` and the beam variants
inside `picket_fence_two_stream_RT_pass` cannot be teased apart from the stats table.

### Grouped shares, and the comparison with the 4-pass form

Rank 0, % of the kernel-time total in each column:

| group | prof_tm (tm+beam) ms | % | prof_0922 (4-pass) ms | % |
|---|---|---|---|---|
| **RT total** | **15428.9** | **81.83** | **15160.8** | **82.00** |
| — column sweep (`picket_fence_*`) | 13567.9 | **71.96** | 13352.0 | **72.22** |
| — conduction (`BuildRadWeights`, `ImplicitRadialUpdate`, `BuildAngularCoeffs`, `NewTimeStep`) | 1427.4 | 7.57 | 1392.6 | 7.53 |
| — apply (`rt_apply` / `par_reduce_clip4`) | 433.6 | 2.30 | 416.2 | 2.25 |
| bvals pack/unpack/seam | 1128.2 | 5.98 | 1119.3 | 6.05 |
| MHD fluxes/update/CT | 995.6 | 5.28 | 915.5 | 4.95 |
| pgen (`SourceFunc`, `HydrostaticEquilibrium`) | 441.6 | 2.34 | 434.9 | 2.35 |
| cs coordinate terms (`Gnomonic*`) | 374.8 | 1.99 | 374.5 | 2.03 |
| EOS c2p | 338.5 | 1.80 | 332.4 | 1.80 |
| other / runtime | 148.1 | 0.78 | 151.0 | 0.82 |
| **kernel-time total** | **18855.7** | | **18488.3** | |
| kernel launches | 70091 | | 66491 | |

**The production form does not change the profile's shape.** RT 81.8 % vs 82.0 %, sweep
72.0 % vs 72.2 % — the tm sweep plus the pseudo-spherical beam cost **+1.6 % in absolute
sweep ms** (13567.9 vs 13352.0) and **+2.0 % in total kernel ms**, and they add
**3600 launches over 300 cycles = 12 per cycle = 6 per RK stage**, which is the second
column pass of the transfer-matrix/adding form. Every non-RT group is within a percent of
the 4-pass arm in absolute ms, as it should be.

Load imbalance is worth recording: rank 1's sweep costs **16308.3 ms against rank 0's
13567.9** (rank-1 kernel total 21696.5 ms, RT 83.8 %). The same imbalance exists in the
4-pass arm but is smaller (rank 1 sweep 13039.9 ms, slightly *below* rank 0), so the
spherical+beam form makes the day/night column asymmetry across the two GPUs worse. The
wall clock is set by the slow rank, so this ~20 % rank spread is now the single largest
cheap win available on this configuration — cheaper than anything inside the sweep.

### Cycles/s

| arm | binary | cpu time, 300 cycles | cycles/s | zone-cycles/s |
|---|---|---|---|---|
| prof_0922 `plain` (4-pass, unprofiled) | a544a761 | 23.042 s | 13.02 | 1.024e7 |
| prof_0922 `prof` (4-pass, profiled) | a544a761 | 23.451 s | 12.79 | 1.006e7 |
| **tm_prof_growth `prof_tm` (tm+beam, profiled)** | 5f348164 snapshot | **26.315 s** | **11.40** | **8.966e6** |

11.40 cycles/s against the 13.0 of the `prof_0922` reference, i.e. **-12.4 %**. Two
things are mixed into that number and should not be: only **+2.0 %** of it is GPU kernel
time (the extra column pass), the rest is host-side — kernel-sum/loop-time falls from
0.79 to **0.717** (18.86 s of kernels inside 26.32 s of loop), which is the extra launch
count and the worse rank imbalance, not arithmetic. And the two rows use **different
binaries**: the snapshot predates the tm-default and FOP-tag commits, worth ~9 %, so at
HEAD the production form should land near **12.4 cycles/s** and the true tm-vs-4-pass
penalty is a few percent, not 12 %. A same-binary A/B (one run with `ck_sweep_form = 0`)
would settle it for ~40 s of apudev; it was not run here.

---

## (8) Split-lag growth gate from a cold start — `g_cfl03` (11941317) / `g_cfl015` (11941318)

The gate the convection boxes used (memory `fmode-dt-taper-test`): an operator-split
radiation lag drives a spurious mode whose growth rate is proportional to the hydro
timestep, so halving `cfl_number` halves it; a physical instability's growth rate does
not move. Both arms are
`inputs/production/deep_hot_jupiter_cs_prod4.athinput` copied **verbatim** (ck_sweep_form
= 1, ck_spherical, ck_beam_sph, max_eta 5e12, pfloor 1e-5, rad_angular false, polytropic
WB with `wb_cache_every = 0`, `dt_min = 1e-3`) with three edits: `cfl_number` 0.3 vs
**0.15**, and `<output1>`/`<output2>` `dt` → 610 s (~0.002 rot) for the fit. **From
scratch**, no restart, 13 min wall each, both finished cleanly (rc = 0, `RUN_WALL` 787 s
/ 785 s).

Reach: `g_cfl03` 9517 cycles to t = 8.153e4 s (0.267 rot); `g_cfl015` 9969 cycles to
t = 5.221e4 s (0.171 rot). Common span **0 → 5.221e4 s = 0.171 rot**; everything below is
measured on it.

### Growth over the common span

| quantity | `g_cfl03` | `g_cfl015` | ratio 0.3/0.15 |
|---|---|---|---|
| 1-KE (radial), end/mid of span | 94.9× | 83.7× | 1.13 |
| KE total, end/mid | 4.29× | 4.11× | 1.04 |
| ME total, end/mid | 5.47× | 5.47× | 1.00 |
| E total, end/start | 0.99800 | 0.99800 | 1.000 |
| **d ln(1-KE)/dt, 2nd half** | **2.205e-4 /s** | **2.087e-4 /s** | **1.057** |
| d ln(KE)/dt, 2nd half | 6.512e-5 /s | 6.295e-5 /s | 1.034 |
| d ln(ME)/dt, 2nd half | 6.449e-5 /s | 6.371e-5 /s | 1.012 |
| 1-KE amplitude, median over 2nd half | — | — | 0.916 |
| total-E drift over the span | -2.015e-3 | -2.018e-3 | 1.001 |

**The growth rate is dt-independent.** Halving the timestep changes d ln(1-KE)/dt by
**5.7 %**, and the total-KE and ME rates by 3.4 % and 1.2 %. A split-lag mode would have
to halve — a factor 2.0 — and it does not come close; the residual few percent is the
ordinary truncation difference between two CFL numbers. The radial-KE e-folding time,
~4.5e3 s = 0.015 rot, is the same physical number in both arms, and the 1-KE *amplitudes*
agree to 8 % at equal times, so it is not merely the rate that matches but the trajectory.
The energy drift (-0.20 % over 0.17 rot) is likewise identical between the arms, i.e. it
is the cold atmosphere radiating as it settles, not a dt-proportional numerical leak.

This is the cold-start spin-up: from a motionless hydrostatic initial state the whole
kinetic and magnetic energy is being built from zero, so a large d ln/dt in the first
fraction of a rotation is expected and carries no information about a split-lag mode by
itself — the *pair* does, and the pair says physical.

### Timestep history

| | `g_cfl03` | `g_cfl015` |
|---|---|---|
| dt at t = 0 | 28.98 s | 14.49 s |
| dt minimum over the run | 5.123 s | 2.667 s |
| dt, mean of last 10 hst rows | 5.286 s | 2.842 s |
| dt trend over the last 30 rows | -1.1e-4 s/s | -8.2e-5 s/s |
| median dt ratio 0.3/0.15 over the common span | **2.021** (last 20 rows: 2.041) | |

dt falls ~5.5× from its t = 0 value as the flow spins up and then **flattens** — the last
30 samples of `g_cfl03` extrapolate to only ~-3 % per further 0.1 rot. `dt_min = 1e-3`
is never approached; there is **no dt collapse in either arm**. The ratio between the arms
holds at 2.02 for the whole span, which is the useful secondary result: dt is on the pure
**hydro CFL** in both, so the `max_eta = 5e12` Ohmic cap (prod4 change 7) is **not
binding from a cold start** either, and `use_rkg_sts = false` is the right setting.

### Floors and C2P

Counters from `dhj.log` (per-interval), normalised to the 786432-cell domain and the
cycles in the interval:

| counter | `g_cfl03` | `g_cfl015` | trend |
|---|---|---|---|
| `eos_dfloor` | 1.49e8 total, ~1.9e-2 cells/cycle | 1.65e8 total, ~2.0e-2 | rises to ~2-3 % of cells by cycle ~1500, then **flat/falling** |
| `eos_tfloor` (200 K) | 3.55e7 total, ~4.2e-3 | 3.00e7 total, ~7.3e-3 | plateaus at 0.4-0.7 % of cells, **no runaway** |
| `eos_efloor` (pfloor 1e-5) | 8.75e6 total, ~4.7e-4 | 5.32e6 total, ~3.0e-3 | 0.03-0.3 % of cells, small but still drifting up at the end in `g_cfl015` |
| `eos_vceil`, `eos_fail`, `c2p_it`, `fofc`, `efloor_de`, `eos_tclamp`, `eos_tset`, `vceil_de` | **all exactly 0** | **all exactly 0** | — |

Three things to read off. (i) **No C2P failure, no FOFC fallback, no velocity ceiling, no
temperature clamp** in ~10^4 cycles on either arm — the inversion is healthy at
`pfloor = 1e-5` barye. (ii) The **tfloor does fire**, at a few tenths of a percent of
cells per cycle, but it *plateaus* and the rate is not systematically larger in the
smaller-dt arm; this is the near-vacuum top of a cold start touching the 200 K guard, not
the runaway the old 50 K study saw (there dt exploded 10 s → 1e9 s in ~1500 cycles; here
dt is flat at 5.3 s after 9500 cycles). (iii) `eos_efloor` is **not** zero, unlike the
`pfloor = 1e-3` claim carried in the prod3 input header; lowering pfloor 100× to 1e-5 did
not eliminate it from a cold start. It is small and it does not feed a dt collapse, but it
is worth one look in the first production dumps to confirm it decays as the atmosphere
settles rather than sitting as a permanent slab.

### Verdict

**dt-independent → not the split lag.** The growth seen from a cold start is physical
spin-up, and there is no numerical kappa-mechanism signature in the prod4 form's
operator-split radiation apply. Together with a clean dt history (no collapse, never near
`dt_min`, pure hydro CFL so `max_eta = 5e12` does not bind), zero C2P failures and zero
FOFC, **prod4 can proceed on the semi-implicit apply from a cold start.**

Two caveats on that clearance. The gate covers **0.17 rotation**; it rules out a mode
with an e-folding time comparable to or shorter than that (which is what the split lag
was), not a slow one — the first production restart at 0.5 rot should be re-fit the same
way. And the `eos_efloor` drift above deserves one confirmation dump. Neither blocks the
launch.

### Raw files

- `prof_tm/prof.log`, `prof.err`, `prof_kk/rank_{0,1}_kernel_{stats,trace}.csv`,
  `rank_{0,1}_domain_stats.csv`, `rank_{0,1}_agent_info.csv`
- `g_cfl03/`, `g_cfl015/`: `run.log`, `run.err`, `dhj.mhd.hst`, `dhj.log`, `bin/`, `rst/`
- each arm keeps its `submit.sh`, its `deep_hot_jupiter.athinput` and a copy of the binary.
