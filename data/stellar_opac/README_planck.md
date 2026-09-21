# Planck-mean opacity tables (for the `<rad_m1>` energy-exchange term)

Same file format and grid as the Rosseland tables here (`ReadOpacityTable`,
`RosselandTable`): header `nT nD lTmin dlT lDmin dlD` = `217 281 2.6 0.025 -14 0.05`,
then `log10 kappa_P [cm^2/g]`, T slowest, second axis `log10 rho` (`table_rho`).

| file | composition | source |
| --- | --- | --- |
| `planck_gs98_x0.7_z0.014.txt` | X=0.700 Y=0.286 Z=0.014, GS98 metals | LANL TOPS/ATOMIC |
| `planck_he_x0.0_z0.02.txt` | X=0 Y=0.98 Z=0.02, GS98 metals | LANL TOPS/ATOMIC |
| `*_ferg+tops.txt` | same | Ferguson et al. 2005 Planck below log T 4.2 (blend 4.0-4.2), TOPS above |

Produced 2026-09-21 with `planck_tools/` (`fetch_tops.py` queries
https://aphysics2.lanl.gov by script; `convert_ferguson.py` reads the Wichita State
`f05.g98.pl.tar.gz` set, interpolating in log Z for Z=0.014; raw downloads are kept in
`bench/m1_opac/`, not in the repository).  TOPS data are valid for log T 3.76-7.07 and
are edge-filled outside.

Checks: the same source's Rosseland mean agrees with the Rosseland tables used so far to
a median 0.005 dex (0.02 dex at the Fe bump); TOPS and Ferguson Planck means agree to a
median 0.00 dex (90 % within 0.32 dex) over log T 3.8-4.5.  Against
`rosseland_gs98_x0.7_z0.014.txt` at log rho = -9: kappa_P/kappa_R = 847 (log T 4.0),
435 (4.7), 38 (5.3), 0.24 (6.0).

Caveats: the Planck mean is ABSORPTION ONLY (no electron scattering; it falls to
1e-5-1e-7 cm^2/g at log T 6.5-7): use it for emission/absorption only, never as a total
extinction.  It is line dominated and non-monotonic at low density; densities below
~1e-8 g/cc may lie under the tabulated minimum for some elements (LANL FAQ).

The tables themselves are fetched/built, not tracked (the `.gitignore` policy of this
directory, see `PROVENANCE.md`); master copies are in `bench/m1_opac/` on viper.  md5:

```
19ff20de3bf3067e7325328e405e6657  planck_gs98_x0.7_z0.014.txt
0cff4fdce01e2d3987d5d87bd12281fc  planck_gs98_x0.7_z0.014_ferg+tops.txt
4cde344b8071184fa8518dd65a88069b  planck_he_x0.0_z0.02.txt
0608d29800a876dfffe13631aa6c5079  planck_he_x0.0_z0.02_ferg+tops.txt
```
