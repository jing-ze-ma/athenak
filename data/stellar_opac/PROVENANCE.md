# Stellar Rosseland opacities (data/stellar_opac)

Third-party opacity tables for the `red_giant` problem generator.  They are fetched, not
tracked (see `.gitignore` here); this file records where each came from and how the
merged table AthenaK reads is built.  All files are for **X = 0.70, Z = 0.02, the
Grevesse & Sauval (1998) mixture**.

## Sources

| file | what | from | md5 |
| --- | --- | --- | --- |
| `oplib_gs98_z0.02_x0.7.data` | Los Alamos OPLIB atomic Rosseland means, OPAL type-1 layout (log T 3.75-9.05, log R -8..1.5) | Farag et al. 2024, Zenodo record 15277019, `share_oplib_type1_tables.tar.gz`, path `kap_data/oplib_gs98_z0.02_x0.7.data` | 1058ab4c321fc6d349c3e17bf73bd87f |
| `aesopus2.0_grains_GS98_Z0.020_X0.7.tab` | AESOPUS 2.0 low-temperature Rosseland means with molecules and solid grains (log T 2.58-4.5, log R -8..1) | Marigo et al. 2022/2023, Zenodo record 8221362, `aesopus2.0_grains_GS98.tgz` | 7696354c251a1f7c0718e6f440a5e2a3 |
| `rosseland_gs98_x0.7_z0.02.txt` | the MERGED table on (log10 T, log10 rho): 217 x 281 nodes, log T 2.6-8.0 step 0.025, log rho -14..0 step 0.05 | `python3 tools/stellar_opac/merge_rosseland.py oplib_gs98_z0.02_x0.7.data aesopus2.0_grains_GS98_Z0.020_X0.7.tab rosseland_gs98_x0.7_z0.02.txt` | (derived) |

Both sources are tabulated against `logR = log10 rho - 3 log10 T + 18`; the merge
evaluates each on the regular (log T, log rho) grid and blends them linearly in log T
between 4.0 (all low-T) and 4.2 (all atomic).  In that band the two agree to a median
0.07 dex, 90th percentile 0.11, max 0.12 -- the atomic table lacks molecules, so below
log T ~ 3.9 it is not to be trusted and is not used.  Nodes outside a source's log R
range take its edge value; a giant envelope lives at log R -3 .. -1.5, inside both.

The LLNL OPAL server (opalopacity.llnl.gov) resets downloads from this cluster; the
Zenodo mirrors above are what actually worked (2026-09-07).

## Why these and not the correlated-k tables

The Exo-FMS correlated-k tables in `data/exo_fms_ck` are premixed exoplanet chemistry on
0.26-325 um bands, valid to 6100 K.  A stellar envelope runs to a million kelvin below
the photosphere and is carried there by the Rosseland mean alone, which is what the
radiative-conduction operator consumes.  Composition: the general EOS should be run with
`eos_xh = 0.70`, `eos_yhe = 0.28` to match.

## Metallicity, and the ck-vs-stellar comparison (2026-09-07)

`rosseland_gs98_x0.7_z0.014.txt` is the same merge at **Z = 0.014**, which is what the
Exo-FMS correlated-k tables are ("1x solar"). The Z = 0.020 table is ~28 % more opaque
along a red giant column. **Use the 0.014 table whenever the run also uses the ck tables.**

**At a red giant photosphere the two datasets AGREE.** Measured along the run's own
column, where the ck table is valid (T 3364-5995 K, tau 1e-4..6, rho 3e-10..2e-8):

    kappa_ck / kappa_stellar(Z=0.014):  median 0.95,  10-90 %  0.84 .. 1.01

Both carry H- bound-free and free-free, which is the dominant continuum there, so this is
what should be expected. No splice is needed and none is kept.
`tools/stellar_opac/merge_ck_stellar.py` remains as a DIAGNOSTIC only.

**Three measurement bugs produced an earlier, wrong claim of a factor 3-5 disagreement.
Anyone repeating this must avoid all three:**

1. the pgen's `column_dump` writes p in **dyn/cm^2, not bar** (a 1e6 error);
2. the ck table's T grid is strongly **NON-UNIFORM** (spacing 0.30 dex at 100 K, 0.015 dex
   at 6100 K), so a lookup must BISECT, as `RosselandTable()` in the code does. Assuming
   uniform spacing was wrong by up to 40x at the hot end;
3. the ck table stops at **6100 K**; deeper points clamp, and the clamped values are not
   a physical disagreement.

Validate any such lookup first: the ck one must reproduce the eight `kappa_R/Freedman`
values the code prints at start-up (it now does, to 1.0000), and the stellar one must
reproduce the pgen's own `kappa_R` column in `column_dump` (it does, to 1.000).

## Where they genuinely differ: GRAINS below ~1500 K

Comparing at the ck table's own nodes (no interpolation, so free of bug 2):

    T < 1500 K    kappa_ck / kappa_stellar = 0.00 - 0.02
    1500-2500 K   0.02 -> 0.36
    T > 2700 K    flat ~0.7 at Z=0.020, ~0.84 at the matched Z=0.014

The low-temperature stellar table is AESOPUS 2.0 **with solid grains**; the correlated-k
tables condense species out of the gas and never restore the dust opacity. This never
touches a red giant (its coldest cell is 3364 K) but it is decisive for brown dwarfs
(photospheres 1000-2000 K) and cold planets.
