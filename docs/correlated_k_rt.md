# Correlated-k radiative transfer (`deep_hot_jupiter_rt`)

Replaces the grey picket-fence two-stream scheme in the `deep_hot_jupiter_rt` problem
generator with **correlated-k**: 11 spectral bands × 8 g-points of premixed equilibrium
opacity, plus CIA, Rayleigh and H⁻ continuum, driven by a blackbody host spectrum.

The motivation is observational. Lee et al. (2021) show that picket fence reproduces
correlated-k *spectra* closely for these atmospheres, but comparing a GCM against
transmission and emission data means the band structure has to be real.

The correlated-k path is **added alongside** the grey one, not in place of it. Runs with
`rt_ck = false` allocate nothing extra, take the original code path, and are bitwise
unchanged — a self-test checks two exact limits of the solver on every correlated-k start,
and the grey path is never touched.

---

## Turning it on

```
<problem>
rt_ck        = true              # false = grey picket fence
ck_pcut_bar  = 10.0              # correlated-k only where p < this
ck_nquad     = 1                 # 1 = diffusivity factor; 2 = 2-point Gauss
ck_star_teff = 6000.0            # host T_eff; the spectrum is a blackbody at this T
ck_table     = data/exo_fms_ck/ck/Premixed_1x_g8_11.txt
ck_data_dir  = data/exo_fms_ck
```

Two things bite immediately:

* **AthenaK command-line overrides only MODIFY parameters that already exist.** A
  `<problem>` block without these lines makes `problem/rt_ck=true` on the command line
  fatal with "not found". `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` lists them all
  with their defaults for exactly this reason.
* **`ck_table` and `ck_data_dir` resolve against the WORKING DIRECTORY**, not against the
  input file. A run launched from a run directory or through `-d` needs absolute paths.

`rt_ck = true` implies `rt_split = true`: the correlated-k solver exists only inside the
split (chain-parallel) kernel. Setting `rt_split = false` alongside it used to run the
grey scheme silently, with every diagnostic looking healthy; it now says what it is doing
and enables the split path.

## The tables are not in git

They come from [Exo-FMS_column_ck](https://github.com/ELeeAstro/Exo-FMS_column_ck) and
carry no upstream licence, so `data/exo_fms_ck/` holds only a `PROVENANCE.md` recording
the download URLs and a `.gitignore`. Fetch them before the first correlated-k run.
`ck_data_dir` must then contain:

```
ck/Premixed_1x_g8_11.txt          the k-table itself (ck_table points at this)
CE_tables/FastChem_ck_1x_int.txt  equilibrium composition: mu, and VMRs of H2 He H e- H-
cia/*_reform_11.txt               four CIA pairs, pre-binned on the 11 bands
ray/Ray_*_11.txt                  Rayleigh cross sections for H2, He, H, e-
sw_flux/sw_band_flux_*_11.txt     stellar spectra (only used when ck_star_teff <= 0)
```

`PROVENANCE.md` also records two ordering traps that were settled empirically, because
the upstream Fortran reader is misleading about both: **the k-table records run DESCENDING
in wavelength while `sw_flux` runs ASCENDING**. The condensation signature in the
0.26–0.42 µm band is the discriminator. A device-side self-test re-checks the band
ordering at every start and aborts if it looks reversed.

---

## What the scheme is

### Bands and g-points

The **Kataria et al. (2013) 11-band grid** — the SPARC/MITgcm binning — with edges at
0.26, 0.42, 0.61, 0.85, 1.32, 2.02, 2.50, 3.50, 4.40, 8.70, 20.00 and 324.68 µm, and 8
g-points per band. Both are fixed by the table (`CK_NB`, `CK_NG`).

The g quadrature is **not** plain Gauss-Legendre. It is split: four nodes weighted to sum
to 0.95 on g ∈ [0, 0.95], then four summing to 0.05 on [0.95, 1]. It must be used exactly
as the table gives it.

This is the same structure Parmentier et al. (2018, A&A 617, A110) use for WASP-121b and
Tan et al. (2024, MNRAS 528, 1016) use across their ultra-hot Jupiter grid, and both keep
the 0.26 µm blue edge even in post-processing.

### The chain count

A "chain" is one column solve. The literature's 88 is bands × g-points; this kernel can
also carry the two-stream angular quadrature, so:

| `ck_nquad` | angular treatment | chains | note |
| --- | --- | --- | --- |
| 1 (default) | diffusivity factor, µ = 1/1.66 | 88 | what GCMs normally do |
| 2 | 2-point Gauss, µ = 0.2113, 0.7887 | 176 | what the picket fence uses |

The approximation in `ck_nquad = 1` is small: after **one** cycle the two differ by 4e-4
in internal energy and 1e-6 in density. Compare schemes at `nlim = 1` — a 200-cycle
comparison measures chaotic divergence, not the scheme.

### Opacity

`k_tot(g, b) = k_ck(g, b) + k_cont(b) + k_Ray(b)`

CIA and Rayleigh are **not** in the k-table but ship pre-binned on the same 11 bands. They
are grey within a band and add to every g-point, so they cost no extra chains. On top of
those, John (1988) H⁻ bound-free and free-free is computed from the FastChem table's
n(H⁻), n(e⁻) and n(H) — it is not optional. Rosseland means of the k-table alone sit at
0.04–0.33 of Freedman et al. (2014) at 3500 K and 0.01–0.06 at 4800 K; with H⁻ and CIA in,
1.04–1.79 and 1.68–4.37. Below 2000 K it does nothing, as expected.

Composition is **premixed** (equilibrium chemistry at fixed metallicity), matched to the
EOS table already in use: [M/H] = 0, H₂ on, ionization on.

### The deep cutoff

Correlated-k runs only where **p < `ck_pcut_bar`** (default 10 bar); deeper than that no
radiative source is applied at all. Two independent reasons, both measured on this setup:

* It is optically thick and convective there. The photosphere (grey τ = 1) sits near
  0.05 bar and the median grey τ is ~9200 at 10 bar, so a Planck bottom boundary at the
  cut is exact to e^-τ.
* The deep interior is 5100–12000 K, far outside any molecular table, and above 5572 K
  more than 1 % of the Planck function falls bluer than the grid's 0.26 µm edge — at
  12000 K it is 30 %. The cut keeps the correlated-k region under ~4800 K, where that
  tail is 3e-3.

The cut is per column (`i_cut(m,k,j)`), since pressure varies day to night. It is where
about a third of the cost saving comes from: p < 10 bar is 65.6 % of cells.

### Boundary conditions

* **Bottom**, at the cutoff: `I_up,b(i_cut) = B_b(T_cut) + I_int,b`, per band.
* **Top**: dark for the longwave; the stellar sweep enters with the cumulative optical
  depth `exp(-tau_down * fac)`.

Band-integrated Planck fractions `f_b(T)` are tabulated on 512 uniform log₁₀T points over
50–20000 K from the Chang & Rhee series, which agrees with direct integration to ~1e-7.
The outer bands are extended to 0 and ∞ so that `sum_b f_b = 1` to 2e-16.

### The stellar spectrum

`ck_star_teff` (default 6000 K) makes the spectrum a **blackbody at that temperature**,
which needs no external data because only the SHAPE is used: the band fractions are
renormalised to sum to 1 and then scaled by the code's own σT_irr⁴, so total insolation
matches the grey scheme and any file's absolute normalisation is irrelevant. Setting
`ck_star_teff <= 0` reads `ck_swflux` from `sw_flux/` instead.

Keep the host in **5500–6500 K**. That is the range the ultra-hot Jupiter literature
works in, and it keeps the flux falling bluer than the 0.26 µm edge — which gets folded
into the bluest band and deposited with near-UV opacities, i.e. too deep — down to 1–3 %.
An A-type host would put 18–32 % outside the grid and would need the 32-band table
instead. Note that `omega` alone does not imply a host: rotation period and T_eq are
independent parameters in this kind of model.

---

## What is deliberately absent

* **A scattering solver.** The two-stream is absorption-only with a scalar albedo; SPARC
  uses Toon et al. (1989). Defensible for clear-sky thermal emission, and it would have
  to change for clouds or hazes.
* **Rayleigh as scattering.** It is added as absorption. Measured contribution to total
  opacity at 3000 K is ≤ 1.4 % worst case and 1e-7…5e-4 in the 0.26–0.85 µm bands that
  carry the stellar flux — metal and molecular lines swamp it at these temperatures. It
  would matter for a cool clear atmosphere.
* **H₂⁻ and He⁻ free-free** (Bell 1980 tables, shipped in `cia/`). Secondary to H⁻, which
  is included. Note the upstream typo in `cia/H2-_ff.txt`: `8.43e02` for `8.43e-2`.

---

## Validation

Four checks, two of them run on every correlated-k start:

1. **Isothermal invariance** (self-test, device, every run). An isothermal atmosphere
   bathed in its own Planck function must carry zero net flux at every level. This works
   because `alp + bet = e0` exactly, so `I = (1-e0) I + e0 B` has B as a fixed point. It
   is the sharpest available check on the recurrence coefficients and the small-x branch.
   Measured: |F|/σT⁴ ≤ 2.6e-17. The run aborts above 1e-12.
2. **Transparent slab** (self-test, device, every run). Optically thin layers over a
   blackbody floor must emit exactly σT⁴. This tests the *weights* rather than the
   recurrence — the g-point weights summing to one, the band Planck fractions summing to
   one, and the flux prefactor being π and not 2π — none of which test 1 can see. Measured
   exactly 1; the run aborts outside 1e-10.
3. **Against Exo-FMS on an identical column.** The longwave agrees to 1.7 %. The
   comparison is made against the production kernel itself through `ck_dump_file`, not
   against a transcription of it.
4. **Global energy.** The base flux is σT_int⁴ to 0.04 % and the shortwave absorbs 98.8 %
   of the analytic incident flux at production resolution.

`tst/test_suite/rad/test_rad_dhj_ck_cpu.py` re-checks 1, 2 and 4 in the regression suite,
and `test_rad_dhj_ck_mpicpu.py` checks that the answer does not depend on the
decomposition. Both build their own binary (the scheme is a user problem generator) and
skip when the tables are absent, which is what CI does.

### Dumping a column

```
ck_dump_file = col.txt           # empty (the default) = off
ck_dump_m    = 0                 # meshblock
ck_dump_j    = -1                # theta index; -1 = centre of the block
ck_dump_k    = -1                # phi index;   -1 = centre of the block
```

Writes level pressures, temperatures, net longwave flux and stellar heating for one radial
column, once, at the first RT call. Note that a value written as `""` in an athinput is
two literal quote characters, not the empty string — leave the value blank to disable it.

---

## Cost

Measured on one MI300A, not extrapolated. The A/B configuration is 8192 columns at
nx1 = 64:

| | integration wall (nlim=500) | all kernels (nlim=100) | RT share |
| --- | --- | --- | --- |
| grey picket fence | 6.434 s | 914.8 ms | 14.5 % |
| correlated-k | 8.623 s = **1.34×** | 1342.1 ms | 40.9 % |

Non-RT work is unchanged (782 vs 794 ms), which is the check that the rewrite left the
rest of the code alone.

**The ratio improves at scale.** Emulating the per-GPU load of a 128 × 64 × 1024 run on 64
nodes — 512 columns per GPU at nx1 = 128 — gives 586 ms grey against 657 ms correlated-k,
i.e. **1.12×**. At that decomposition the *grey* RT is the starved one: 512 threads is 8
wavefronts on 304 CUs, and it is 41.7 % of GPU time on its own. Correlated-k's split
kernel exposes 512 × 22 threads instead and the extra 84 chains ride on idle hardware.
Real multi-node runs add communication to both schemes, diluting the ratio further toward
1, so 1.12× is an upper bound. The extra arrays are ~40 MB/GPU there.

### Notes for anyone optimising it further

`rt_chain` is ~95 % of the RT cost and is bound by the recurrence dependency chain,
`I_down[i] <- I_down[i+1]` through an `expm1`, with the whole grid resident at once — so
there are no spare waves to hide latency behind. The following have all been **measured
and rejected**: flipping the k-table layout, blocking the lookup by band, caching layer
coefficients between sweeps, dropping the private intensity column, `RT_NB` of 1, 2 or 8,
precomputing kappa per (cell, chain) (tried on both array layouts), and skipping `expm1`
where the layer is optically thick.

The one change that did work was memory coalescing, and it was worth 3.4×. `par_for`
flattens its **last** argument fastest, so in
`par_for(..., 0,nmb1, 0,nblk-1, ks,ke, js,je, lambda(m,blk,k,j))` adjacent threads differ
in `j`; the per-cell arrays had to be indexed `(m, slot, i, k, j)` rather than
`(m, slot, k, j, i)` to put adjacent lanes in the same cache line. **Check the
thread-to-index mapping against the array layout before believing any null result on this
kernel** — several of the rejections above first looked like proof that loads were not the
limit, and were re-measured after the fix.

`RT_FP32` (compile-time, default off) makes the recurrence single precision: 1.42× on
`rt_chain`, but that is 25 % of RT and only 6 % of the run, and it costs 1.7e-4 on the net
longwave flux. Judge this kernel on the total, not on the RT.

The radial array size is chosen at run time from a tiered set (72, 136, 264, 520), so
`nx1` up to 512 needs no recompile. `nx1 = 512` itself still fails, but on a pre-existing
core limit — Kokkos cannot find a valid team size for the MHD flux kernels, whose LDS
scratch scales with `ncells1` — not on anything in the RT.
