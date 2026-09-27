---
name: dhj-overall-gpu-profile
description: "Where GPU time actually goes in a deep-hot-Jupiter run (job 10980944): mhd::MHD::CalculateFluxes is 50.6%, RT 14.5% grey / 41% correlated-k, and everything else in the pgen is under 8%."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T22:30:00.000Z
---

## THE ONE TO USE: xe_long-matched production profile (2026-08-22, MI300A)

Matched to **`/viper/u2/jinma/ATHENAK/bench/xe_long/b10_e13`**, the 2500 K hot Jupiter long
runs: 64 x 64 x 128 in 32 meshblocks, x1 = 9.44e9..13.04e9, nghost 2, stretch 2.0,
**bbot = 10 G, max_eta = 1e13, NO RKG**, pfloor 1.0e0, dfloor 7.26e-12,
`ohmic_resistivity = eos` -- but with `eos = general/table` and `rt_ck = true`.
Reproduce from the shipped `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` with

    mesh/x1max=13.04e9 meshblock/nx2=16 meshblock/nx3=16 mhd/max_eta=1.0e13 \
    mhd/use_rkg_sts=false mhd/pfloor=1.0e0 mhd/dfloor=7.26e-12 problem/bbot=1.0e1 \
    problem/rt_ck=true problem/ck_table=... problem/ck_data_dir=...

20 cycles, 339 ms of GPU kernel time:

| share | kernel group |
|---|---|
| **37.6 %** | correlated-k RT (chain kernel **30.7 %**, rt_pre + rt_apply 7.0 %) |
| **32.5 %** | MHD CalculateFluxes |
| 6.3 % | boundary buffers / halo / prolongation |
| 5.9 % | other pgen source terms |
| **5.4 %** | ConsToPrim -- the tabulated EOS inversion |
| 3.3 % | RK update / PrimToCons / corner-E / dt |
| 3.1 % | resistivity (eos x_e, no STS) |
| 2.3 % | user hydrostatic BC |
| 3.6 % | unclassified |

### max_eta = 1e13 vs 1e14, same long-run setting (2026-08-22, one job, one node)

200 cycles FROM THE INITIAL CONDITION, tabulated EOS + correlated-k:

| max_eta | RKG STS | cpu (s) | t reached | dt |
|---|---|---|---|---|
| 1e13 | off | 4.385 | 2682.445 | 13.41 |
| 1e14 | off | 4.361 | 2682.445 | 13.41 |
| 1e14 | **on** | **5.692** | 2682.445 | 13.41 |

**At t = 0 the cap makes NO difference at all** -- identical dt, identical cost. dt_diff is
1.93e4 against a hydro dt of 13.41, so tau ~ 7e-4 and diffusion is nowhere near binding
with the initial bbot = 10 G field. **Turning RKG on there costs a flat 1.30x for nothing**,
which is the "pure ~40 % loss when diffusion is not binding" the input file warns about,
measured. The substage count is floored at 3.

**RESOLVED by running it out.** Three 4-minute runs, SEQUENTIAL on the same MI300A
(never co-scheduled -- two jobs on one node share the host and bias the timing), ~10,000
cycles each, t ~ 3.2-3.7e4:

| config | t reached | cycles | mean dt | sim-s per wall-s | vs 1e13 |
|---|---|---|---|---|---|
| **1e13, no RKG** | 3.721e4 | 10670 | 3.487 | **158.8** | **1.00x** |
| 1e14, no RKG | 3.243e4 | 10650 | 3.045 | 138.3 | 0.87x |
| 1e14, RKG on | 3.170e4 | 8176 | 3.877 | 135.3 | 0.85x |

dt against time (t = 5e3, 1e4, 2e4, 3e4, final):

| config | | | | | |
|---|---|---|---|---|---|
| 1e13, no RKG | 13.41 | 4.45 | 2.88 | 2.29 | 2.17 |
| 1e14, no RKG | 13.41 | 4.45 | 2.89 | **1.582** | **1.582** |
| 1e14, RKG on | 13.41 | 4.45 | 2.88 | 2.32 | 2.25 |

All three track each other until ~t = 2e4, then **1e14 LOCKS at dt = 1.582** -- the diffusive
limit binding. That is the SAME 1.582031 [[xe-resistivity-long-runs]] recorded for the
ideal-gas runs, which makes sense: it depends only on max_eta and the grid, not on the EOS
or the RT.

**RKG does exactly what it claims and is still not worth it.** It frees dt back to 2.25
(1.42x over the locked 1.582) but costs 1.30x per cycle, so it lands at 0.85x -- marginally
WORSE than just leaving it off. **Use max_eta = 1e13 with RKG off.**

Caveat: t ~ 3.2e4 is just past the crossover. The ideal-gas runs at t ~ 4.5e4 found 1e13
ahead by 2.5x at 10 G, so the gap keeps widening; 1.15x is a lower bound, not the asymptote.
Extrapolating from the final dt alone (2.17 vs 1.58) gives 1.37x at equal cost per cycle.

Kernel shares at 1e14 with RKG on shift only because the extra substages add cheap work:
RT 30.9 %, fluxes 26.7 %, ConsToPrim 10.0 %, boundary 12.1 %, resistivity 4.6 %.
Without RKG, 1e13 and 1e14 are identical to 0.1 %.

### Cost of the new physics against what the long runs actually ran

200 cycles from the IC, same grid and floors, current binary:

| configuration | zone-cyc/s | dt | cpu-s per sim-s | vs the long runs |
|---|---|---|---|---|
| ideal + grey (what xe_long ran) | 4.029e7 | 10.14 | 1.283e-3 | 1.00x |
| tabulated EOS + grey RT | 3.116e7 | 13.41 | 1.255e-3 | **0.98x** |
| tabulated EOS + correlated-k | 2.395e7 | 13.41 | 1.632e-3 | **1.27x** |

**The tabulated EOS is FREE here** -- 1.29x more per cell exactly cancelled by a 1.32x larger
timestep. All of the 1.27x is the correlated-k RT (1.30x on its own).

Caveat: dt = 10-13 is the value from the INITIAL CONDITION. `xe_long` recorded mean dt 3.19
at t ~ 1e5 once the flow developed, so absolute cpu-s per sim-s will rise for every row; the
RATIOS should hold. Also note the xe_long table's 2.49e7 zone-cycles/s is NOT comparable --
that binary predates the RT coalescing fix and this session's EOS work.

## Earlier: same profile at nx1=128, max_eta=1e14, RKG on, bbot=5 (NOT the long-run setup)

`rocprofv3 --kernel-trace --stats`, one MI300A (apudev vipa1001), 128 x 64 x 128 in 32
meshblocks (8192 columns, nx1=128), **tabulated general EOS + correlated-k RT +
ohmic_resistivity=eos**, max_eta=1e14 and RKG STS ON (the input's own production values),
20 cycles, 764 ms of GPU kernel time.

| share | kernel group |
|---|---|
| **35.7 %** | correlated-k RT (rt_chain_ck alone **30.2 %**; rt_pre + rt_apply 5.5 %) |
| **27.8 %** | MHD CalculateFluxes |
| **12.2 %** | ConsToPrim -- the tabulated EOS inversion |
| 9.6 % | boundary buffers, halo exchange, prolongation |
| 4.0 % | other pgen source terms |
| 3.5 % | resistivity + RKG super-time-stepping |
| 2.5 % | user hydrostatic BC |
| 2.1 % | RK update / PrimToCons / corner-E / dt |
| 2.6 % | unclassified |

**RT + fluxes are 63.5 % between them.** That is where any further optimisation has to go.
The EOS is 12.2 %, which is the ceiling on all remaining EOS work -- see
[[general-eos-optimization]].

This is KERNEL time, not wall: host gaps and launch overhead are excluded, and there are
1879 launches in the boundary group alone, so its wall share is larger than 9.6 %.

Raw CSVs: `/viper/ptmp/jinma/claude_eos_gpu/prof/prod_kernel_stats.csv`, job script
`prof.sh` alongside. Kernel names are mangled Kokkos templates -- match on the enclosing
function (`picket_fence_two_stream_RT`, `CalculateFluxes`, `GeneralMHD::ConsToPrim`), not
on the par_for label, which is NOT in the symbol.

---

## Earlier profile (grey RT, kept for the grey comparison)

Full kernel trace, correctness.athinput (meshblock 64x16x16, 32 blocks, 1 APU), nlim=100,
ms per 100 cycles. Grey picket fence in the left column, correlated-k on the right.

| kernel | grey | share | correlated-k |
|---|---|---|---|
| **`mhd::MHD::CalculateFluxes`** (par_for_outer) | **463.4** | **50.6 %** | 463.8 |
| `picket_fence_two_stream_RT` (pgen) | 132.6 | 14.5 % | 553.6 |
| `SourceFunc` / `usrsource` (pgen) | 42.6 | 4.7 % | 45.9 |
| `mhd::MHD::RKUpdate` | 36.0 | 3.9 % | 36.0 |
| `HydrostaticEquilibrium` (pgen user BC) | 28.4 | 3.1 % | 28.8 |
| `IdealMHD::ConsToPrim` | 19.2 | 2.1 % | 19.3 |
| `Coordinates::SrcTermsSphericalPolarMHD` | 18.0 | 2.0 % | 18.0 |
| `Resistivity::AddEMFGeneralResist` + `AddFluxGeneralResist` | 26.5 | 2.9 % | 26.5 |
| boundary pack/unpack CC + FC | ~49 | 5.4 % | ~52 |
| TOTAL GPU | 914.9 | | 1346.2 |

**The single biggest kernel in the whole run is the MHD flux calculation, not the RT** --
3.5x the RT in the grey configuration, and still 34 % of GPU time with correlated-k on. If
the run as a whole needs to be faster, that is where the time is.

**Nothing in the pgen outside the RT is significant.** `usrsource` (42.6 ms) is a per-cell
EOS pressure plus TempKelvin, which are table lookups on the general EOS -- that is why it
costs what it does. It is par_for over (m,k,j,i) with i fastest, matching the (m,n,k,j,i)
array layout, so it is properly coalesced; no pathology. The user BC is 28.4 ms over 603
launches. Together 7.8 %.

## CalculateFluxes, looked at (job 10980959, 10980975)

It is three `par_for_outer` team kernels in `src/mhd/mhd_fluxes.cpp`, HLLD
(`MHD_RSolver` 3), split by direction:

| kernel | teams | LDS | ms/100cyc |
|---|---|---|---|
| `mhd_flux1` (x1), teams over (m,k,j) | 10368 | 10824 B | 137.2 |
| `mhd_flux2` (x2), teams over (m,k), j serial inside | 576 | 15200 B | 166.0 |
| `mhd_flux3` (x3), teams over (m,j), k serial inside | 576 | 15200 B | 160.1 |

Counters: x1 occupancy 10.37 waves/CU, VALUBusy 13.2 %; x2/x3 occupancy 7.63, VALUBusy
8.9 %. **VALUUtilization is only 53 % in all three** -- half the lanes in every wave idle --
and the GPU is ~90 % idle in the kernel that is half the runtime.

**x2 and x3 have 18x fewer teams than x1** because they are launched over two indices and
sweep the third serially inside the team, and they cost 326 ms against x1's 137. That is
where the headroom is, and it is the same shape of problem as the RT had.

**Tested and REJECTED: it is not wavefront alignment.** `ncells1 = nx1 + 2*nghost` is 68,
not a multiple of 64, which looked like the obvious cause of 53 %. Per-cell flux cost:
nx1=60 (ncells1=64) 7.520, nx1=64 (68) 7.214, nx1=124 (128) 5.820. Aligning to exactly 64
is WORSE. What the trend actually shows is a fixed per-team cost amortised over longer
i-ranges -- and that is not actionable here, since the meshblock already spans the whole
x1 extent.

**The real lever turned out to be the TEAM SIZE, not the x2/x3 sweep** -- see
[[par-for-outer-team-size]]. 1.24x on the flux kernels and 1.11x whole-run at this grid,
bitwise identical, but it is core code and was deliberately NOT applied.

Note the meshblock decomposition is already tuned -- see [[meshblock-decomposition-gpu]].
