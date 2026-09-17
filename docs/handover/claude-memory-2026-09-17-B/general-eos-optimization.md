---
name: general-eos-optimization
description: "MEASURED 2026-08-22: the tabulated general EOS was 62% of a dhj run; five changes (dcab8ff7, 8c75d780) made it 1.41x faster overall. Next big idea is a direct T(rho,e) inverse table, best evaluated on GPU."
metadata:
  node_type: memory
  type: project
---

## How to measure it at all

`general_eos=gamma` keeps the ENTIRE general code path and swaps only the table, so
`table` vs `gamma` isolates the table cost exactly. Both are run-time knobs -- no rebuild.
Compare **zone-cycles/cpu_second**, not wall, since the modes take different dt. Add
`mhd/ohmic_resistivity=perna` or `gamma` mode fatals ("found no electron fraction").

Baseline on dhj 64x8x8, nlim=100, serial CPU: **table 7.218e4, gamma 1.883e5 -> the table
was 62 % of the run.**

An `eos_dlog` scan 0.025 -> 0.2 (61 MB -> 0.4 MB table) moves throughput **0.7 %**, so it
is NOT footprint-bound. That one scan is the cheapest diagnostic here; run it first.

## What was done (7.218e4 -> 1.016e5, 1.41x)

| change | commit | gain |
|---|---|---|
| log-domain Newton residual | dcab8ff7 | **1.28x** |
| skip the unused mu interpolation in PressureAndGamma1 | dcab8ff7 | +1.1 % |
| Hermite derivatives templated on demand; SolveLog/EvalResidual templated on mode | 8c75d780 | +1.9 % |
| fused TemperaturePressureGamma1 (solve and evaluate share their logs) | 8c75d780 | **+4.8 %** |
| tbl transposed (ITNVAR,ny,nx) -> (ny,nx,ITNVAR) | 8c75d780 | +1.9 % |

**The big one:** the surfaces are log10 of the SPECIFIC quantity, so
`log10(e) = log10(rho) + E(x,y)` identically and the residual is a subtraction. The old
code exponentiated z to get T, logged it again inside Interpolate, exponentiated the
interpolated value, and logged that -- 2 exp + 3 log10 per iteration computing the
identity, plus log10(rho) recomputed for a fixed density.

Correctness: bitwise on dhj (radiation off AND on) and on 3 of the 4 shipped
`linear_wave_*_geneos*`; `linear_wave_mhd_geneos_table` moves 0.02-0.12 % in L1 with the
convergence rate 1.965-1.968 against its 1.8 threshold.

## General EOS vs IDEAL, on the CURRENT deep_hot_jupiter_rt (measured 2026-08-22)

`deep_hot_jupiter_rt_ideal_xe.athinput` vs `..._eos.athinput`, same pgen, same grid, and
the settings that are NOT the EOS forced to match:
`mhd/max_eta=1.0e13 mhd/use_rkg_sts=false problem/bbot=5.0` (the two inputs ship with
different values for all three). Both keep `ohmic_resistivity = eos`, so the x_e table
lookup is common to both and only the thermodynamics differ. Grey RT on in both.

| | per zone-cycle | per unit of SIMULATED TIME |
|---|---|---|
| before this session (9a8fe0bd) | **2.81x** | 2.08x |
| after dcab8ff7 + 8c75d780 | **2.01x** | **1.49x** |

ideal 2.717e5, general 1.354e5 zone-cycles/s. **The general EOS also takes a 1.35x LARGER
timestep** -- lower Gamma_1 in the ionization zones means a lower sound speed -- so per
unit of physical time the penalty is 1.49x, not 2.01x. Always quote both; zone-cycles/s
alone overstates the cost by a third.

That 2x is the WHOLE application (MHD + grey RT + resistivity), not the EOS kernel.

## GPU RE-MEASUREMENT (2026-08-22, MI300A, apudev vipa1001) -- READ THIS BEFORE OPTIMISING FURTHER

**The CPU numbers above overstate the EOS's importance by about 5x.** Same A/B, one MI300A,
grey RT, 128 x 64 x 128 in 32 meshblocks (8192 columns), nlim=200. Run-to-run noise 0.25-0.6 %.

| | CPU (4096 cells) | GPU (1.05M cells) |
|---|---|---|
| table share of the run, before | **62 %** | **18.1 %** |
| table share of the run, after | 46 % | **12.5 %** |
| table evaluation itself | 1.89x faster | **1.54x faster** |
| WHOLE RUN | **1.41x** | **1.062x** |

So the optimisations are real and transfer (the table really is ~1.5-1.9x cheaper), but on
the GPU the EOS was never the bottleneck -- `CalculateFluxes` and the RT are. **The GPU gain
is 6 %, not 41 %.** With correlated-k on it is 1.5 %.

**Prediction that FAILED:** I expected the node-major transpose to pay MORE on GPU than on
CPU (24 cache lines -> 4 per Eval). The table speedup is 1.89x on CPU and only 1.54x on GPU,
so it paid LESS. Do not assume a coalescing argument transfers; measure it.

**CEILING ON FURTHER EOS WORK: 12.5 % of the run** (15.6 % at the 64-node per-GPU load).
A direct T(rho,e) inverse table might take two thirds of that, i.e. ~8 %. **Not a priority.**

Caveat: on GPU the `table` runs used the input's own `ohmic_resistivity=eos` while the
`gamma` control used `perna`, so the quoted share includes the x_e lookup as well as the
thermodynamics -- the thermodynamics alone is even smaller.

## GENERAL vs IDEAL, ON GPU -- the number to quote

| per-GPU load | per zone-cycle | dt | **per unit of SIMULATED TIME** |
|---|---|---|---|
| 8192 columns (nx1=128) | 1.45x | 1.32x larger | **1.10x** |
| 512 columns (64-node production) | 1.38x | 1.32x larger | **1.04x** |

**On the GPU the general EOS costs about 10 % of wall time, and only 4 % at production node
counts** -- against 1.49x on CPU. It gets CHEAPER relatively as the per-GPU load drops.
The 1.32x larger timestep (lower Gamma_1 through the ionization zones -> lower sound speed)
is doing most of the work; never quote zone-cycles/s alone.

Also measured: correlated-k costs **1.36x** of grey per simulated second at 8192 columns,
confirming the 1.34x in [[correlated-k-design]].

**How to run these:** binaries and job scripts in `/viper/ptmp/jinma/claude_eos_gpu`
(`gpu.sh`, `gpu2.sh`). Two traps: the session scratch under `/tmp` is NODE-LOCAL, so stage
everything on `/viper/ptmp`; and **`export HSA_XNACK=1` is REQUIRED** for the
`Kokkos_ARCH_AMD_GFX942_APU` build or every run dies with "Memory access fault by GPU"
after Kokkos warns about xnack. `apudev` has 2 nodes, 15 min limit.

## Facts worth keeping

- **2.94 residual evaluations per root find.** The warm start works as documented; there
  is nothing to win in the iteration count on CPU.
- **~7.7 root finds per cell per cycle** (2 RK stages + the RKG STS substages + ghost
  work). This multiplier, not the per-solve cost, is what makes the EOS expensive. NOT
  yet investigated -- see below.
- `mu` is wanted only by resistivity and the cooling source term, never by ConsToPrim.

## Next, in order

1. **A direct T(rho,e) inverse surface.** e is monotonic in T at fixed rho, so tabulate
   log10 T on a (log rho, log e) grid at build time and replace the 2.94-patch Newton
   solve with ONE interpolation, or with one interpolation plus a single Newton step for
   exactness. **On GPU this is worth much more than on CPU**, because a wavefront pays
   max(iterations) not mean -- the variable trip count is pure divergence today.
2. **Cut the 7.7 c2p per cell-cycle.** Check whether every RKG super-time-stepping
   substage really needs a full temperature solve.
3. **fp32 table storage.** Halves the traffic; a node becomes one 64 B cache line.
4. Anything cache-related must be judged on a GPU. See [[correlated-k-design]]'s coalescing
   lesson: check the thread-to-index mapping against the array layout before believing any
   null result.
