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
