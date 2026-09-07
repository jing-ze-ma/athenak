# Exo-FMS correlated-k data, 11-band (Kataria+2013) grid

Downloaded 2026-08-21 from <https://github.com/ELeeAstro/Exo-FMS_column_ck> (`main`),
the code companion to:

> Lee, E. K. H., Parmentier, V., Hammond, M., et al. 2021, MNRAS 506, 2695,
> "Simulating gas giant exoplanet atmospheres with Exo-FMS: comparing semigrey,
> picket fence, and correlated-k radiative-transfer schemes", arXiv:2106.11664

**The upstream repository carries no LICENSE file.** These files are redistributed here
for research use only. Cite Lee et al. (2021) and the underlying line lists, and check
with the author before publishing results that depend on them.

## Contents

| path | what |
| --- | --- |
| `wavelengths_GCM_11.txt` | 12 band edges in um, DESCENDING: 324.68 ... 0.26 |
| `ck/Premixed_1x_g8_11.txt` | premixed k-table, 1x solar, equilibrium condensation |
| `cia/*_reform_11.txt` | H2-H2, H2-He, H2-H, He-H collision-induced absorption |
| `cia/H2-_ff.txt`, `cia/He-_ff.txt` | H2- and He- free-free |
| `ray/Ray_*_11.txt` | Rayleigh cross sections, H2 / He / H / e- |
| `sw_flux/sw_band_flux_{W121,HD189}_11.txt` | stellar flux per band, EXAMPLES ONLY |

## k-table format (HELIOS-k / `ck_form == 2` in `src/ck_opacity_mod.f90`)

```
line 1                     species list (text)
                           nT nP nband ng        -> 38 34 11 8
                           T(1..nT)              [K]    100 .. 6100
                           P(1..nP)              [bar]  1e-8 .. 1000
                           wl(1..nband+1)        [um]   descending, 324.68 .. 0.26
                           wn(1..nband+1)        [cm-1] ascending
                           Gx(1..ng)  Gy(1..ng)  g nodes and weights
                           kappa: for iT, for iP, for b = nband..1 : ng values
```

Everything after line 1 is whitespace-separated numbers, so it can be read with a single
token stream.

**Ordering, verified against the data (2026-08-21).** The k-table records run in the same
order as the `wl` edges, i.e. DESCENDING wavelength: the first record of each (T,P) block
is 324.68-20 um and the last is 0.26-0.42 um.

This is the opposite of what Exo-FMS's own reader appears to do (`do b = nwl, 1, -1` in
`ck_opacity_mod.f90`), so do not take that loop at face value. The discriminator is
condensation. With band 11 = 0.26-0.42 um, band-mean kappa at 0.1 bar is

| T [K] | 300 | 800 | 1500 | 2500 | 3500 |
|---|---|---|---|---|---|
| 0.26-0.42 um | 1.3e-8 | 1.3e-7 | 0.71 | 31 | 60 |
| 20-324.68 um | 14 | 14 | 3.5 | 3.8 | 0.64 |

The optical cliff between 800 and 2500 K is TiO/VO/Fe/Na/K coming out of condensation; the
far-IR decline is the H2O rotational band. Reversed, both are physically impossible.

**`sw_flux` is ordered ASCENDING in wavelength, opposite to the k-table and to `wl`.**
Reversing the W121 file and normalising reproduces a 6460 K blackbody to about 1 % per
band (0.2283 vs 0.2260 at 0.61-0.85 um, 0.1906 vs 0.1905 at 0.85-1.32 um); as listed it is
exactly backwards.

**Units:** kappa is cgs, cm^2/g. Exo-FMS multiplies by 0.1 only to convert to MKS; AthenaK
is cgs throughout, so the values are used as read. Interpolate log10(kappa) in
(log10 T, log10 p).

**Continuum:** CIA and Rayleigh are NOT in the k-table. They are grey within a band and
add to every g-point: `k_tot(g,b) = k_ck(g,b) + k_cont(b) + k_Ray(b)`.

## Known defects in the upstream files (checked 2026-08-23)

The files here are byte-for-byte upstream. Two things in them are wrong, and are recorded
rather than silently patched, so that this directory stays a faithful copy of what was
downloaded.

1. **`cia/H2-_ff.txt`, row `4142`, second column reads `8.43e02` where it should be
   `8.43e-2`** — a factor of 1e4, obvious from its neighbours `5.84e-2` and `1.01e-1` and
   from the smooth Theta dependence of every other row. **It does not affect any run:**
   the continuum reader opens only `cia/{H2-H2,H2-He,H2-H,He-H}_reform_11.txt`,
   `ray/Ray_{H2,He,H,e-}_11.txt`, `CE_tables/FastChem_ck_1x_int.txt` and the k-table.
   `H2-_ff.txt` and `He-_ff.txt` ship unused — H2- and He- free-free are deliberately not
   modelled (see docs/correlated_k_rt.md). Fix the value if that ever changes.

2. **`wavelengths_GCM_11.txt` is a count line `12` followed by ASCENDING edges**
   (0.260 ... 324.68), not the descending list the table above describes. Also unused: the
   band edges are read from the k-table itself.

Everything else was checked and is clean:

| check | result |
| --- | --- |
| k-table value count | 113696 = 38 x 34 x 11 x 8, exact |
| k-table monotonic in g (required of a k-distribution) | all 14212 chains pass |
| k-table negatives | none; range 2.6e-175 .. 2.0e4 |
| FastChem records | 8194 as declared, no leftover tokens |
| FastChem mu | 0.8006 .. 2.3267, no negative VMRs, none summing > 1.5 |
| `cia/He-_ff.txt`, CIA and Rayleigh tables | no order-of-magnitude outliers beyond the 1e-99 no-data sentinels and genuine Rayleigh falloff |
