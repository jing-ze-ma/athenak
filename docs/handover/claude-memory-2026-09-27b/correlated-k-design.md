---
name: correlated-k-design
description: "KERNEL IS WIRED (42690945) and correlated-k RUNS: FULL SCHEME (longwave + shortwave) RUNS at 2196 ms/100cyc, wall ~2.4x (5x worse than the 409 ms projected: 15.5 ms/chain not 5.4). Physics validation ladder NOT yet done. Replacing the dhj picket-fence RT with correlated-k. Decisions settled 2026-08-21: premixed 1x solar, Kataria 11-band, Lee Exo-FMS tables, 10 bar cutoff. Tables downloaded to data/exo_fms_ck/ and the reader is COMMITTED (52e29b12). Next: band-integrated Planck fractions, then the kernel."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T19:00:00.000Z
---

## RESUME HERE (2026-08-22)

**The correlated-k RT is BUILT, VALIDATED and OPTIMISED.** Branch `polar-average-perf` at
**35c8483b**, pushed to the fork. Turn it on with `problem/rt_ck=true`; it implies `rt_split` and picks
its own radial array size. Tables live in `data/exo_fms_ck/` (fetched, not in git --
`PROVENANCE.md` says where from).

State: 11 Kataria bands x 8 g-points x 1 angular point = 88 column solves, diffusivity
factor 1/1.66, correlated-k above 10 bar only, blackbody host at 6000 K, CIA + Rayleigh +
H- continuum, premixed 1x solar FastChem chemistry. **557 ms/100cyc RT, 1.34x total wall
vs the grey picket fence at the A/B config and ~1.12x at 64-node per-GPU load.**

Validated: longwave within 1.7 % of Exo-FMS on an identical column, base flux =
sigma T_int^4 to 0.04 %, shortwave absorbs 98.8 % of analytic, two exact-limit self-tests
pass at machine precision every run.

**NOTHING IS LEFT.** Physics done, optimised, documented (`docs/correlated_k_rt.md`),
regression-tested (`_cpu` + `_mpicpu`), MPI-verified, harness retired. The only open
question was ever the host T_eff, and that is settled at 5500-6500 K.

Note `bench/polar_ab/fixed200.a` is STALE (predates the b4e0953c stellar-heating fix); the
live reference is `fixed200_swfix.a` and all A/B scripts point at it.

---

Cost side is settled in [[rt-chain-parallel-split]]; this is the physics/implementation
side.

## Decisions (all the user's, 2026-08-21)

1. **Science driver: observations.** The user heard from Vivien Parmentier that
   correlated-k is needed for observational comparison. This overrides Lee+2021's finding
   that picket fence reproduces c-k spectra closely. Do not relitigate it.
2. **Premixed tables** (equilibrium chemistry, fixed metallicity) -- not on-the-fly
   mixing. Match the composition to the EOS table already in use: [M/H] = 0, H2 on,
   ionization on.
3. **Deep cutoff at 10 bar**: correlated-k only where p < 10 bar, something grey
   (Rosseland diffusion) below.
4. **Band structure: Kataria+2013 11-band grid** (the SPARC/MITgcm binning, 0.26-324.86
   um). **Tables: Lee's Exo-FMS tables.**

## Tables: FOUND, DOWNLOADED, VERIFIED

From <https://github.com/ELeeAstro/Exo-FMS_column_ck> into `data/exo_fms_ck/` (NOT tracked
in git -- no upstream licence; `PROVENANCE.md` and a `.gitignore` are tracked instead, and
they record the download URLs). The band grid IS Kataria's, so both decisions hold.

`ck/Premixed_1x_g8_11.txt`: 38 T x 34 p x 11 bands x 8 g. **T 100-6100 K, p 1e-8 to
1000 bar**, 1x solar, equilibrium condensation, 33 species incl. TiO VO Fe FeII Na K.
10x/100x/1000x and 32-band variants also exist upstream. Band edges are exactly Kataria:
0.26 0.42 0.61 0.85 1.32 2.02 2.50 3.50 4.40 8.70 20.00 324.68 um.

The T ceiling of 6100 K covers the 4822 K the 10 bar cutoff implies, so **the 10 bar
cutoff stands** -- no need to retreat to 3 bar.

**The g quadrature is NOT plain Gauss-Legendre.** Split: 4 nodes weighted to sum 0.95 on
g in [0,0.95], then 4 summing to 0.05 on [0.95,1]. Must be used exactly as given.

**CIA and Rayleigh are NOT in the k-table** but ARE provided pre-binned on the 11-band grid
(`cia/`, `ray/`). They are grey within a band and add to EVERY g-point:
`k_tot(g,b) = k_ck(g,b) + k_cont(b) + k_Ray(b)`. Not extra chains.

**Two ordering traps, both settled empirically -- see PROVENANCE.md for the evidence:**
k-table records run DESCENDING in wavelength (opposite to what Exo-FMS's own reader
suggests), while `sw_flux` runs ASCENDING. Do not trust the Fortran loop; the condensation
signature in the 0.26-0.42 um band is the discriminator.

**Still missing:** a stellar flux file for the actual host star. Only W121 and HD189 ship
with the repo. Bin a PHOENIX/Kurucz spectrum onto the same 12 edges; also pin down the
normalisation (the W121 values sum to ~9.8e7 cgs, which is not obviously sigma T_irr^4).

## Measured facts this rests on (from `fixed200.a`, 524288 cells)

- Photosphere (grey tau = 1) is at **~0.05 bar**. Median tau: 3.3 at 0.1 bar, 166 at
  1 bar, **9180 at 10 bar**. The cutoff is ~4 decades below the photosphere, so a Planck
  bottom boundary there is exact to e^-tau.
- p < 10 bar is **65.6 %** of cells -> RT ~409 ms/100cyc (from 623 for the full column),
  wall ~1.2x. Cutoff trade: 0.3 bar 46.9 %/292 ms/3714 K, 1 bar 53.0 %/330 ms/3752 K,
  3 bar 57.8 %/360 ms/3853 K, 10 bar 65.6 %/409 ms/4822 K.
- Whole domain is 4e-6 to 238 bar, T median 3737 K, 90th percentile 10113 K. **40.6 % of
  all cells exceed 4000 K**, which is where `get_kapr` (Freedman+2014) clamps -- the deep
  region is 5100-12000 K and already extrapolated today. The cutoff isolates that problem
  rather than solving it, which is fine: below 10 bar it is optically thick and convective.

## COST REALITY CHECK (job 10979675) -- the projection was 8x optimistic

Correlated-k longwave, 176 chains, 10 bar cut: **3367 ms/100cyc** (rt_pre 640, rt_chain
2712, rt_apply 15) against the grey picket fence's 115. Span goes ~1530 -> ~4780 ms, so
**wall ~3.1x**, not the 1.2x quoted during design. Two compounding errors, both in the
earlier estimate:

1. **176 chains, not 88.** The literature's "88 column solves" is band x g-point. This
   kernel ALSO carries the 2-point Gauss angular quadrature the picket fence uses
   (mug = 0.2113, 0.7887). Factor 2.
2. **15.5 ms/chain, not 5.4.** The synthetic harness modelled the table lookup but not the
   extra `exp` undoing log10 kappa, the per-band Planck source read, or the continuum
   read. Factor 2.9.

The 10 bar cut delivered its 0.66 as predicted.

**Optimisations tried and their real value:** removing `pow()` from H- and restricting the
continuum to above the cut took rt_pre 1042 -> 640 ms (bitwise unchanged). **Flipping the
k-table layout to (b,g,iT,iP) so the (T,p) plane is contiguous did NOTHING** -- 2725 ->
2711 ms, 0.5 %. The k-table access pattern is not the bottleneck, which is the same answer
the synthetic blocked-lookup experiment gave. Do not try it a third time.

**DONE: `problem/ck_nquad` = 1 (default) uses the diffusivity factor mu = 1/1.66; 2 keeps
the 2-point Gauss.** 176 -> 88 chains, RT 3260 -> **2045 ms/100cyc** (rt_pre 640,
rt_chain 1396, rt_apply 8), **wall ~2.3x**. The approximation is small: after ONE cycle
the two differ by 4e-4 in eint and 1e-6 in density. Always compare schemes at nlim=1 --
a 200-cycle comparison measures chaotic divergence, not the scheme.

## SPEED, after optimisation (25bc700a, b018a9dc)

**RT: 2234 -> 557 ms/100cyc, 4.0x.** Two changes: rt_pre parallelised (1.34x) and the array
layout fixed for coalescing (3.4x on rt_chain).

**TOTAL COST, measured (job 10980774), not extrapolated -- earlier wall-time figures in
this thread were extrapolations from a remembered span and were WRONG:**

| | integration wall (nlim=500) | all kernels (nlim=100) | RT share |
|---|---|---|---|
| grey picket fence | 6.434 s | 914.8 ms | 14.5 % |
| correlated-k FP64 | **8.623 s = 1.34x** | 1342.1 ms | 40.9 % |
| correlated-k FP32 | 8.104 s = 1.26x | 1209.8 ms | 34.4 % |

Non-RT work is unchanged, 782 vs 794 ms -- the check that the rewrite left the rest alone.

**Part one, rt_pre:** It was doing per-CELL work
(continuum, band Planck, table index) in a kernel parallel over (m,k,j) only: 0.52 waves/CU
and 2.8 % VALUBusy, the same starvation the monolithic RT had. Split into four kernels --
mu0 per column, T/p per cell, cut index per column, opacity/Planck per cell -- and it went
**635 -> 68 ms**. Verified by column dump: bitwise identical output.

**rt_chain is now 95 % of the RT at 1588 ms and has resisted six attacks, ALL MEASURED:**

| tried | result |
|---|---|
| k-table layout (T,p) plane contiguous per chain | 0.5 % |
| blocking the lookup by band, sharing (T,p) weights | 0.4 % |
| caching layer coefficients between sweeps | **-120 %** |
| dropping the private intensity column (scratch) | **-22 %** |
| more chain blocks, RT_NB 2 / 1 | **-16 % / -24 %** |
| precomputing kappa per (cell,chain) in a parallel kernel | -1 % |
| FP32 recurrence (`RT_FP32`, kept, default off) | +2.8 %, costs 1.7e-4 on the flux |

The last rules out almost everything: removing the four table loads, the exp AND the
interpolation from the inner loop did not move it. Counters say 10.85 waves/CU, VALUBusy
4.85 %, MemUnitStalled 1.65 %, and at 704 workgroups **the whole grid is resident at once**,
so there are no spare waves to hide latency. It is bound by the recurrence dependency
chain, I_down[i] <- I_down[i+1] through an expm1.

**Two more tried on the fixed layout, both NULL (do not retry):**
- Skipping `expm1` where the layer is optically thick (`x > 40`, where e0 is exactly 1):
  458 vs 452 ms, slightly WORSE. Too few chain-cells qualify, or the divergence costs more.
- Precomputing kappa per (cell,chain), RETESTED now that its extra load is coalesced:
  533 vs 452 ms, 18 % WORSE. The 613 MB array's traffic beats halving the kappa work.
  This one has now failed on BOTH layouts.

**The GREY path is already at its best**: rt_split=false 132.7 ms vs rt_split=true 167.4.
With only 4 chains the split path gains no parallelism and pays extra launches and global
flux traffic. Do not switch production to the split path.

**FP32 (`RT_FP32`, default OFF).** On the STARVED layout it bought 2.8 %; on the fixed
layout it buys **1.42x on rt_chain**, 457 -> 325 ms, i.e. 25 % of RT. Still off, because
25 % of RT is only **6 % of the run** (8.62 -> 8.10 s) and it costs 1.7e-4 on the net
longwave flux. Judge this kernel on the TOTAL, not on RT.

## IT WAS UNCOALESCED MEMORY. 3.4x (626f50f4)

`par_for` flattens its LAST argument fastest. In
`par_for("rt_chain_ck", ..., 0,nmb1, 0,nblk-1, ks,ke, js,je, lambda(m,blk,k,j))`
adjacent threads differ in **j**. The per-cell arrays were `(m,slot,k,j,i)` with i fastest,
putting adjacent LANES n1 = 68 doubles = **544 bytes** apart: every lane its own cache
line, 8 useful bytes per 64 fetched.

Re-indexing them **`(m,slot,i,k,j)`** so j is fastest: **rt_chain 1586 -> 461 ms**, RT total
1670 -> **557 ms/100cyc**, column dump BITWISE IDENTICAL. rt_apply 14 -> 24 and rt_pre_opac
64 -> 67, noise against 1125 ms saved.

**Several of the seven verdicts were artefacts of the starvation, re-measured after the
fix (b42a9a62):** FP32 +2.8 % -> +1.42x; RT_NB=2 -16 % -> faster on rt_chain but slower on
the total (nblk doubles, so rt_apply sums twice as many slots -- RT_NB stays 4); dropping
the private column -22 % -> neutral. Genuinely useless on both layouts: k-table layout,
blocking the lookup by band, precomputing kappa. **Do not trust a null result on this
kernel without first checking memory is not in the way.**

**Why the seven earlier attempts all failed, and why one of them misled me:** they were all
aimed at compute, occupancy or scratch. The kappa-precompute null result looked like proof
that loads were not the limit -- it was not: `ck_lk` is 888 kB and one (band,g) plane is
10 kB, so those four table loads were CACHE HITS. The expensive loads were `Bb_g` and
`kc_g`, scattered across 76 MB. Cache residency, not load count. **Lesson: on this kernel,
check the thread-to-index mapping against the array layout before anything else.**

**Old note:** rt_pre is ~635 of the 2196, nearly a third, and it is per-cell continuum and
Planck setup rather than column solves. It does NOT scale with the chain count. Note it ALSO still runs the
grey optical-depth sweep and the old 3-band Q_v when correlated-k is on -- dead work.

## Build list, in order of pain

**DONE so far. Nothing is wired into the RT kernel yet; `rt_ck=false` is bitwise
unchanged and that is checked every run.**

- **52e29b12** table ingestion. `problem/rt_ck`, `problem/ck_table`,
  `problem/ck_pcut_bar`. `read_ck_table()` fills `ck_lk(iT,iP,band,g)` (log10 kappa, g
  fastest), `ck_lT`, `ck_lP`, `ck_gw`, `ck_wl`, validating monotonic grids and weights.
- **2d4d45cb** band-integrated Planck fractions. `build_planck_fractions()` tabulates
  f_b(T) on 512 uniform log10 T points, 50-20000 K, from the Chang & Rhee series (agrees
  with direct integration to ~1e-7). The 0.26-324.68 um grid does not capture the whole
  Planck function, so the outer bands are extended to 0 and infinity and sum_b f_b = 1 to
  2e-16. **Above 5572 K more than 1 % of the flux is bluer than the grid** -- the 10 bar
  cut keeps this setup under 4822 K where the tail is 3e-3; the deep interior at 12000 K
  would be 30 %, so the BAND STRUCTURE, not just the table range, needs that cut.
  `ck_planck_bands()` is the device-side lookup.
- **fb713a74** `ck_tp_index()` (clamped bracketing search -- clamp the FRACTION too, that
  was the old NaN bug) and `ck_kappa()` (bilinear in log10 kappa). Plus `ck_selftest()`,
  which runs on device at startup and fatals if the band ordering looks reversed, using
  the condensation signature. It reproduces an independent Python parse digit for digit.

- **(this commit)** continuum. `read_ck_continuum()` reads the FastChem CE table
  (`CE_tables/FastChem_ck_1x_int.txt`, SIX columns per record: mu then VMR of H2 He H e-
  H-), the four CIA pair tables and the Rayleigh cross sections. `ck_continuum()` returns
  kappa_cont for all 11 bands at a cell. Matches an independent parse to every digit.
  Needs `problem/ck_data_dir` as an ABSOLUTE path -- runs happen in subdirectories.

- **(H- commit)** John (1988) H- bound-free and free-free, using n(H-), n(e-), n(H) from
  the FastChem table. Rosseland/Freedman ratios go from 0.04-0.33 to 1.04-1.79 at 3500 K
  and 0.01-0.06 to 1.68-4.37 at 4800 K. Below 2000 K it does nothing, as expected.
  **The opacity side is now complete enough to drive the kernel.**

**Validation done offline:** table Rosseland means match Freedman+2014 to a factor
1.1-2.7 over 1000-2500 K. At 3500 K they diverge because a Rosseland mean is set by the
transparent windows and **CIA and H- free-free are not in the k-table** -- that is the
next piece, and it is not optional.

0. **H2- and He- free-free** from the Bell (1980) tables in `cia/`. Secondary to H-,
   which is now in. NOTE an upstream typo in `cia/H2-_ff.txt`: `8.43e02` for `8.43e-2`.
1. **DONE. The shortwave.** Folded into the longwave down-sweep (they share the kappa
   lookup, so a separate kernel would pay twice): full scheme 2196 ms/100cyc vs 2045
   longwave-only, so the shortwave costs 7 %. `problem/ck_swflux` picks the stellar
   spectrum, default WASP-121 -- **REPLACE with the real host star**. Only the SHAPE is
   used; values are renormalised to sum to 1 and scaled by the code's sigma T_irr^4, so
   total insolation matches the grey scheme and the file's absolute normalisation never
   matters. The file is ordered ascending in wavelength, opposite to everything else.
   Energy closes: most transparent (band,g) reaches tau = 162 at the 10 bar cut, surviving
   flux 1.5e-72 of incident. The longwave needs nothing structural:
   `dtau_i = tau[i]-tau[i+1]` is exactly `kapr*rho*dr`, purely local, so per-chain kappa is
   a lookup and no new recurrence. The stellar sweep uses CUMULATIVE tau
   (`exp(-gamv*tau_down_r_f[i]*fac)`), so each stellar (band,g) gains its own downward
   recurrence and must move out of `rt_pre` into a chain-parallel kernel. Needs the host
   star's spectrum binned on the same bands.
2. **Band-integrated Planck fractions** B_b(T), 11 per cell instead of sigma T^4/pi,
   tabulated in T. Lands in `rt_pre`, currently the 67 ms floor -- expect it to double.
3. **Bottom BC at the cutoff**: `I_up,b(i_cut) = B_b(T_cut) + I_int,b`, replacing
   `I_up[is] = Iint + I_down[is]`. `Iint` becomes per-band B_b(T_int).
4. **Per-column cutoff index**: pressure varies day to night, so `rt_pre` computes
   `i_cut(m,k,j)` and `rt_chain` sweeps `[i_cut, ie+1]`. Neighbouring columns are similar,
   so intra-wave imbalance is mild.
5. **Scattering**: the current two-stream is absorption-only with a scalar albedo. SPARC
   uses Toon+1989 with scattering. Defensible to skip for clear-sky thermal emission;
   Rayleigh in the visible is not negligible and clouds would force it.
6. **Validation ladder**: 1 band x 1 g-point with a grey kappa must reproduce today's
   answer; then a line-by-line 1D column comparison on a fixed T-p profile; then global
   energy balance against sigma T_eff^4.

## SCALING: correlated-k gets CHEAPER, relatively, at high node count (job 10981278)

Emulating the per-GPU load of a 128 x 64 x 1024 run on 64 Viper nodes (128 GPUs), i.e.
**512 columns/GPU with nx1 = 128**, against the 8192-column config everything else was
measured on:

| per-GPU load | grey RT | c-k RT | all GPU grey | all GPU c-k | ratio |
|---|---|---|---|---|---|
| 8192 columns, nx1=64 | 132.7 | 554.9 | 917.3 | 1347.7 | **1.47x** |
| **512 columns, nx1=128** | 244.4 | 307.9 | 586.3 | 656.8 | **1.12x** |

**At high node count correlated-k costs only ~12 % more GPU time, not 47 %.** With 512
columns the GREY RT is catastrophically starved -- 512 threads is 8 wavefronts on 304 CUs,
and its time even goes UP (133 -> 244 ms) despite 1/16 the columns, because a fully
latency-bound kernel is priced by per-thread work and the column doubled. Correlated-k's
split kernel exposes 512 x 22 = 11264 threads and the extra 84 chains ride on idle
hardware. Grey RT is then 41.7 % of GPU time on its own -- "RT is cheap" fails for the
picket fence too at that scale.

Caveat: single-GPU emulation of the per-GPU load, no MPI. Real 64-node runs add
communication to BOTH schemes, which dilutes the ratio further toward 1. So 1.12x is an
upper bound. Extra c-k arrays are ~40 MB/GPU at that decomposition.

**The grey path cannot be rescued the same way**: with only 4 chains there is no chain
dimension to parallelise over.

## Radial size is now automatic (d05acd08)

RT_NNC used to be a fixed 72 and correlated-k FATALLED at nx1 = 128. Sizing it generously
is not free -- chain kernel at n1=68 costs 454 ms at NNC 72, 503 at 136, 540 at 272. So the
kernel is now instantiated at **72, 136, 264, 520** and the smallest fitting n1 is
dispatched at run time via a generic lambda taking `std::integral_constant`. hipcc accepts
a device lambda inside one. Bitwise unchanged, no measurable cost, nx1=128 just works.
RT_NNC still guards the GREY split path only.

**The same was done for the MONOLITHIC grey RT, which is the PRODUCTION path.** It had a
flat `constexpr int NN = 270` with NO check, while the long runs use nx1 = 256 -> n1 = 260:
ten cells from silently overrunning five per-thread arrays, and nghost = 4 would have gone
over. Now tiered 72/136/264/520/1032. Bitwise identical, cost unchanged.

**nx1 = 512 still does not run**, but that is a PRE-EXISTING core limit, verified on the
pre-dispatch binary: Kokkos cannot find a valid team size for the MHD flux kernels, whose
LDS scratch scales with ncells1. The RT is no longer what binds first.

The flat 270 also survives in `double_gray_two_stream_RT` and its `_source` variant -- both
dead, never called.

## WHAT IS LEFT (reviewed 2026-08-22)

**Blocking real science:**
1. **OPEN DECISION: which host star, and therefore which band grid.** `problem/ck_star_teff`
   (7a390ea2) now builds the spectrum as a blackbody at a given host T_eff -- no external
   data, since only the SHAPE is used. What is missing is the T_eff itself.

   **The setup implies an A star.** No host is specified anywhere (the scheme only ever
   needed T_eq), but omega = 2.06e-5 is a 3.53 day period, and if tidally locked,
   T_irr = T_eff sqrt(R/a) gives T_eff = 9200 K (2 Msun) to 12300 K (0.8 Msun).

   **That breaks the 11-band choice.** Kataria's grid stops at 0.26 um, fine for
   HD 209458 at 6100 K (1.3 % of flux bluer) but NOT for an A star:

   | host T_eff | flux bluer than 0.26 um |
   |---|---|
   | 5800 | 1.3 % |
   | 6446 (the W121 placeholder) | 2.6 % |
   | 10000 | **18.3 %** |
   | 12300 | **31.8 %** |

   That flux is folded into the bluest band and deposited with near-UV opacities, so it
   lands too deep. Alternatives shipped in the same repo: **32-band, blue edge 0.20 um,
   6.7 % outside at 10000 K, premixed table AVAILABLE** (`Premixed_1x_g8_32.txt`), which
   is 32 x 8 = 256 chains, about 2.9x the RT cost of 88; and 34-band, 0.115 um, ~0 %
   outside, but NO premixed table shipped.

   **RESOLVED by the literature, see [[uhj-band-structure-literature]]:** Parmentier+2018
   (WASP-121b) and Tan+2024 both use exactly 11 Kataria bins x 8 k-coefficients for
   ultra-hot Jupiters, and both keep the 0.26 um blue edge even in post-processing. All
   their hosts are 5500-6500 K. And omega = 2.06e-5 does NOT force an A star -- Tan+2024
   vary T_eq and rotation period as INDEPENDENT parameters, which is standard. So: keep the
   11-band grid, set `ck_star_teff` in 5500-6500 K, and the blue tail is 1-3 %.

**Quantified and judged FINE for this planet -- do not spend time on these:**
2. **Rayleigh is added as ABSORPTION, not scattering** (there is no scattering solver at
   all). Measured contribution to total opacity at 3000 K: <= 1.4 % worst case, and
   1e-7..5e-4 in the 0.26-0.85 um bands that carry the stellar flux. Metal and molecular
   lines swamp it at these temperatures. Would matter for a COOL clear atmosphere, or with
   clouds.
3. **H2- and He- free-free** (Bell 1980 tables in `cia/`) are not included. Secondary to
   H-, which is in. NB `cia/H2-_ff.txt` has an upstream typo, `8.43e02` for `8.43e-2`.
4. Sub-cycling: ruled out on physics, see the t_rad/dt measurement.

**Deliberately absent:** a scattering solver. Needed only for clouds/hazes; SPARC uses
Toon+1989. Absorption-only is defensible for clear-sky thermal emission here.

**Process, per CLAUDE.md / CONTRIBUTING:**
5. **DONE (5cb3e1bb).** `tst/test_suite/rad/test_rad_dhj_ck_cpu.py`. Passes in 73 s on
   viper13 (~60 s of that is the build at -j256). Two awkwardnesses, both unavoidable:
   the scheme is a USER pgen so the test compiles its OWN binary into `tst/build_ck`
   (the driver's shared build cannot carry a `PROBLEM=`), and it `pytest.skip`s when
   `data/exo_fms_ck/` is missing, which is what CI will do. Checks: the two device
   self-tests, sum_b f_b = 1, "88 column solves", that rt_ck took the split path, and two
   physical ones on a 64x8x8 grid -- nightside F(icut)/sigma T_int^4 = 1.068 (tol 20 %)
   and dayside absorbed/incident = 0.965 (tol 10 %); production gives 1.0004 and 0.988.
   Nightside column is `ck_dump_k=2` (mu0 = -0.918), dayside `ck_dump_k=5` (+0.918).
   The same commit adds `ck_dump_file`/`ck_dump_m`/`ck_dump_j`/`ck_dump_k` to the input
   file so the overrides resolve. **`ck_dump_file = ""` in an athinput sets the name to
   two literal quote characters** -- `ParseLine` does not strip quotes -- so the default
   must be written as a bare `ck_dump_file =` (affa34e6 fixed that; the first version
   dumped a column into a file called `""` on every run).

   A second test, `test_rad_dhj_ck_mpicpu.py`, checks decomposition independence; shared
   build/run/parse helpers live in `test_suite/rad/dhj_ck_common.py`, which is not a test
   module. Each test compiles its OWN binary (`tst/build_ck`, `tst/build_ck_mpi`) and
   removes it in a `finally`; ~85 s each.

   **A plain CPU serial build WORKS again** (g++ 11.5, no modules, ~60 s at -j256) --
   `bfield_bcs.cpp`'s MPI_Allreduce is now guarded. [[viper-hip-build-recipe]] said
   otherwise and was stale.

   **Style note:** `python run_test_suite.py --style` currently FAILS on this checkout for
   reasons unrelated to any of this -- cpplint hits on `src/bvals/*` and flake8 hits on
   untracked `run/*.py` and `tools/solar_convection/*.py`. Do not read that as your fault.
6. **DONE (1261459e): `docs/correlated_k_rt.md`**, in-repo alongside `general_eos.md`
   rather than on the GitHub wiki. Covers the scheme, the options, the table layout and
   its two ordering traps, what is deliberately absent, the validation ladder, the cost
   table, and the full list of rejected optimisations with the warning about null results
   on an uncoalesced layout. The pgen used to point at `bench/exofms_compare/README.md`,
   which is NOT in the tree (`bench/` was never committed -- so `bench/polar_ab/*.a`
   referenced elsewhere in this note is scratch, not a repo path); it now points here.
7. **DONE 2026-08-22 (affa34e6). MPI is CLEAN, but getting there found a real bug.**
   c-k on 1/2/4/8 ranks (8 meshblocks, 10 cycles): `.hst` AND `.bin` bitwise identical to
   serial. Grey too. The per-rank k-table read is a non-issue: setup is 5.6 s at every
   rank count and the GREY run, which reads no tables, takes the same -- the 15 MB
   general-EOS table build every rank does anyway swamps the 3.1 MB of c-k tables.

   **The bug was NOT in the RT: `mhd_corner_e.cpp` / `_uct.cpp` left the polar EMF
   average's host mirrors (`polar_inner_h`, `polar_outer_h`) at extent 0**, because
   006bae07 guarded their sizing on `inner_local.extent(0) != ncells1` and mhd.cpp's ctor
   already pre-sizes `inner_local`. Only the MPI branch touches the mirrors, so serial
   never noticed. **Every multi-rank run with `use_polar_boundary` died on cycle 1**
   ("deep_copy extents of views don't match: (0) inner(68)") from 2026-08-19 to the fix
   in **3882e37f**. Loud, so no wrong answers -- but no multi-rank dhj run worked in that
   window. There is now an `_mpicpu` test (below) and it fails on the parent commit.

   Test harness note: `mpirun` needs `module load gcc/14 openmpi/5.0` -- MPCDF modules are
   hierarchical, so there is no `mpicxx` in a bare shell. `pytest`/`flake8` need
   `module load python-waterboa/2025.06`; the system python3 has neither.
8. **DONE (35c8483b). `rt_nchain` and `rt_ktab` are GONE**; `rt_nchain` survives only as
   a derived count (4 grey, CK_NB*CK_NG*ck_nquad c-k). That removed the per-chain weight
   multiply, the `synth` flag and the dead (T,p) proxies from both inner sweeps, a 422 kB
   synthetic table and a weight vector allocated every run, and the "RT scaling harness
   active ... table lookup off" line that was flatly wrong under rt_ck.

   **The `rt_split` A/B was KEPT** -- it is the check that the chain-parallel kernel
   reproduces the monolithic one, and it was not even reachable from the shipped input
   (overrides only modify existing parameters), so `rt_split` is now listed there too.
   Verified bitwise on all three paths -- grey monolithic, grey split, correlated-k --
   history output and the c-k column dump, plus both regression tests.

**Also: the correlated-k options are now listed in `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`.**
Not cosmetic -- AthenaK overrides only MODIFY existing parameters, so without them
`problem/rt_ck=true` fatals with "not found". `ck_table`/`ck_data_dir` resolve against the
WORKING DIRECTORY, so runs from a run directory need absolute paths.

**Fixed 2026-08-22:** `rt_ck=true` with `rt_split=false` used to silently run the GREY
scheme with every diagnostic looking healthy. rt_ck now implies rt_split.

## Kernel config to build into

Split path (`problem/rt_split`), **RT_NB = 4 to 8**, RT_CACHE off. Cost is ~5.4 ms/100cyc
per chain including the table lookup, and blocking the lookup by band buys nothing, so
(nband, ng) is a free choice on physics grounds. `RT_NNC` must be >= n1.
