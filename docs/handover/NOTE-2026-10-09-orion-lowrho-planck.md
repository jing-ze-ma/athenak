# NOTE 2026-10-09 (orion): LTE Planck and Rosseland means at very low density (AG Car thin atmosphere)

Task: `TASK-2026-10-09-orion-lowrho-planck.md` (rt-integration 0f6a02ed). CPU only, no AthenaK source change.
Deliverables are in `docs/handover/agcar-opac-orion-1009/` (`tables/`, `plots/`, `scripts/`). The work dir is
`/orion/ptmp/jinma/agcar_opac_1009` (raw downloads, per-T npz, logs). SLURM job 213657 (1 orion node, 233 runs,
about 1 h). No restart was used because this is a table calculation.

## Result in one paragraph
The LTE Planck mean at rho 1e-20..1e-15 g/cm^3 is now computed, not extrapolated. Relative to viper's ext2
extension it is within +-0.1 dex at 2800-3500 K and **up to +0.65 dex higher** at 3800-5600 K (Fe II, Mg II,
Ca II stay important). At 5600-7000 K and rho <= 1e-18 it is up to **-0.85 dex lower** (Mg/Fe already doubly
ionised in LTE). At 7500-8000 K and rho <= 1e-19 it is **+0.8..+0.9 dex higher**. At 10-20 kK it is
**-0.2..-0.5 dex lower**, and down to **-1.0 dex** at 20 kK and rho <= 1e-19, where ext2 carries the flat
TOPS clamp. The
Rosseland mean is electron scattering plus a small line excess, and it agrees with ext2 to <= 0.12 dex
everywhere (median 0.00). Validation: the Rosseland mean matches Ferguson 2005 and TOPS to 0.1 dex (max).
The Planck mean matches Ferguson to 0.06 dex at T <= 3500 K and TOPS to 0.06-0.15 dex at 6300-7900 K. It is **0.3-0.55 dex below
BOTH references at log T >= 4.05** (11-20 kK). That deficit is NOT explained by the line list, the H level
cutoff, broadening or the frequency grid (all tested, section 4). Treat kP above 11 kK as a lower bound by
~0.5 dex. The 1-2 dex uncertainty of the ext2 Planck extension therefore comes down to about +-0.3 dex below
10 kK and about -0/+0.5 dex above.

## 1. Method (`scripts/`)
- **Mixture.** GS98 metal number fractions (the LANL `grevsau1` list, the same as TOPS and
  `data/stellar_opac/planck_tools/mixspec.py`), X 0.36, Y 0.62, Z 0.02. There are 26 elements (H, He, C-Zn without
  Li-B and F; F is 3e-5 of the metals by number).
- **Ionisation.** Own Saha chain for every element, all stages up to XI. The ionisation energies come from NIST
  ASD. The partition functions are sums over NIST ASD levels below the ionisation limit. H I and He II are
  hydrogenic with n <= 30, and n_max 15 or 60 changes nothing (section 4). H- is in Saha equilibrium. There is
  no continuum lowering, which is irrelevant at n_e <= 1e12. Against FastChem 4.0.3 (gas, with ions) at
  2800-5000 K, n_e agrees to <= 0.02 dex.
- **Molecules.** pyfastchem 4.0.3 (the RSG venv `/orion/ptmp/jinma/rsg_wind_1008/grains/venv`, JANAF logK, no
  condensation) for log T <= 3.95. The atoms bound in CO, CN, OH, SiO, H2O, TiO and H2 are removed before Saha.
  Molecular line opacity uses the DACE 1e-8 bar per-species cross sections (0.01 cm^-1, the RSG slab
  `/orion/ptmp/jinma/rsg_wind_1008/dace/raw`, read in place). The slab stops at 2900 K, so it is held at 2900 K
  above that; OH goes to 4900 K. CO dominates: at 2818 K and 1e-12 essentially all C is in CO.
- **Lines (bound-bound).** The lists are Kurucz CD-ROM 1 (`lowlines` + `highlines`, 42 M lines including the
  predicted Fe-group lines) for every ion where it has more lines, and Kurucz `gfall08oct17` for the rest. That
  covers H I, the light-element ions and Sc-Ni ions with few lines (`merge_lines.py`). Isotope and hfs
  fractions are applied. The line strength is S = (pi e^2/m_e c) gf exp(-E_l/kT)(1-exp(-h nu/kT))/U.
- **Planck mean** = sum over lines of S B_nu(nu_l)/B (exact; it needs no frequency grid and no profile), plus
  continuum, plus molecules integrated on their native 0.01 cm^-1 grid. It is absorption only.
- **Bound-free.** The ground configuration uses the Verner et al. 1996 analytic fits (`phfit2.f`, shifted by the
  level energy). All other levels use hydrogenic Kramers with Z_eff = ion charge + 1, n_eff from the binding
  energy and g_bf = 1 (an "OP-like hydrogenic" treatment without resonances). H I uses Verner for n=1 and
  Kramers for n >= 2; He II is hydrogenic. Free-free is hydrogenic (g_ff = 1) for all ions. H- bf and ff use
  the Gray (2005) fits. Rayleigh scattering (H, He, H2, Kurucz/Dalgarno fits) and Thomson scattering are in
  the Rosseland mean only. CIA, H2- and He- ff are omitted because they scale as rho^2 or n_e n_He, which is
  negligible below 1e-12.
- **Rosseland mean.** Opacity sampling on a log-uniform grid from 300 to 4e5 cm^-1 (33 um-25 nm; it captures
  >= 99.98 % of the B and dB/dT weight at every T) at R = nu/dnu = 1e6. The profiles are Voigt (Humlicek W4),
  with Doppler width from thermal + 2 km/s microturbulence (the Ferguson convention) and radiative damping
  (gfall Gamma_R, classical if missing). The wings are cut at 1 % of nu. Lines with peak opacity < 1e-10
  cm^2/g at their maximum population are skipped in the grid only; they remain in the exact Planck sum.

## 2. Broadening statement
At rho <= 1e-12 (n_H <= 4e11, n_e <= 3e11 cm^-3), pressure broadening is negligible. Van der Waals gives
Gamma_6 ~ 10^-7.5 n_H <~ 1e4 s^-1, and quadratic Stark gives Gamma_4 ~ 1e-5 n_e <~ 3e6 s^-1. Both are well
below the radiative Gamma_R ~ 1e7-1e9 s^-1, and Gamma_R itself is below the Doppler width (a = Gamma/4 pi
dnu_D ~ 1e-3-1e-2). The lines are Doppler cores with radiative-damping wings. Pressure broadening is not
included.
- **DACE.** The 1e-8 bar molecular cross sections are Doppler-limited: the Lorentz HWHM is ~1e-9 cm^-1 against a
  Doppler HWHM of ~0.04 cm^-1. Our gas at 2800 K and 1e-12 sits at ~1.6e-7 bar, also Doppler-limited. DACE's
  lowest-P point is therefore NOT pressure-broadened beyond Doppler.
- **Profiles were built from the Kurucz list for atoms**, not taken from DACE.
- **Sensitivity.** The Planck mean is independent of broadening by construction (exact line sum). The Rosseland
  mean changes by <= 0.011 dex with no microturbulence and by <= 0.004 dex with the wing cut at 0.3 % or 3 %.

## 3. Frequency grid and convergence (`tables/convergence.txt`, `plots/convergence.png`)
The sampled Planck mean (opacity sampling, the same grid as the Rosseland mean) is compared with R = 1e6:

| R | kP sampled, max / median abs dlog [dex] | kR, max abs dlog [dex] |
|---|---|---|
| 1e4 | 1.04 / 0.25 | 0.002 |
| 3e4 | 0.29 / 0.08 | 0.0004 |
| 1e5 | 0.0018 / 0.0000 | < 1e-4 |
| 3e5 | 0.0006 / 0.0000 | < 1e-4 |

At R = 1e6 the sampled Planck mean equals the exact line sum to 0.0005 dex at every grid point. The Doppler
width is b/c ~ 7e-6 for Fe at 2800 K with 2 km/s, so sampling needs R >~ 1e5. The Rosseland mean is
continuum (Thomson) dominated and converged at R = 1e4. The production tables use the exact sum for kP and
R = 1e6 for kR.

## 4. Validation (`tables/validation_*.txt`, `tables/summary.txt`, `plots/*_cuts.png`)
The overlap is Ferguson 2005 GS98 (Wichita `f05.g98.pl` / `f05.gs98`, X 0.35/0.5 interpolated to 0.36; the same
files as viper's `sources/ferguson05`) at log R >= -8, and TOPS (x0.36 z0.02, the viper table) at log T >= 3.764
and log rho >= -14.

| quantity | vs Ferguson 2005 (N=1689) | vs TOPS (N=902) |
|---|---|---|
| kR | median 0.000, p10/p90 -0.046/+0.012, max 0.080 dex | median -0.013, p10/p90 -0.031/+0.034, max 0.097 |
| kP | median +0.025, p10/p90 -0.37/+0.30, max 0.48 dex | median -0.28, p10/p90 -0.53/+0.06, max 0.85 |
| (reference spread) | Ferguson - TOPS on the same cells: kP median -0.09, p10/p90 -0.36/+0.40; kR median -0.016 | |

kP by temperature (this work minus reference, median over the overlap rows):

| log T | 3.45-3.55 | 3.60-3.85 | 3.90 | 3.95-4.00 | 4.05-4.30 |
|---|---|---|---|---|---|
| vs Ferguson | -0.06..+0.05 | +0.13..+0.34 | +0.01 | -0.26..-0.35 | -0.30..-0.45 |
| vs TOPS | - | -0.06..-0.15 (3.80-3.85) | -0.11 | +0.07..+0.10 (range -0.46..+0.59) | -0.30..-0.55 |

- **3.45-3.55 (2800-3500 K):** agreement with Ferguson to 0.06 dex. Ca II, Mg II and Fe I/II lines dominate;
  CO adds up to 18 %.
- **3.6-3.85:** 0.1-0.35 dex above Ferguson, but within 0.06-0.15 dex of TOPS where TOPS exists (6300-7100 K).
  Ferguson is itself 0.2-0.3 dex below TOPS there, so this is within the reference spread.
- **3.95-4.0 (9000-10000 K):** H recombination. kP varies by > 1 dex per 0.1 dex in T and the references
  disagree by up to 1 dex with each other (TOPS 0.49 vs Ferguson 1.43 at 8900 K and 1e-14; we give 1.08).
  Differences there are T-interpolation dominated.
- **>= 4.05 (11-20 kK): 0.3-0.55 dex below both references. Not explained.** Tested and excluded:
  - predicted Fe-group lines (CD-ROM 1 vs gfall08oct17 only: <= 0.016 dex anywhere);
  - Kurucz's newer computed C III list (`atoms/0602/gf0602z.lines`, 38826 lines vs 8673): +8 % on C III, the
    dominant contributor;
  - H I/He II n_max 15-60: 0.000 dex;
  - bound-free + ff + H-: <= 0.25 % of kP anywhere;
  - line skip threshold 1e-12: 0.000 dex;
  - ionisation: Saha with NIST U. At 20 kK and 1e-12 the dominant ions are C III, Si IV, O III, N III and
    Fe IV, whose strong lines are in the far UV (x = h nu/kT > 7).

  Remaining suspects (not tested): line data of the dominant light ions beyond Kurucz/NIST (OP/ATOMIC
  compute all configurations), photoionisation resonances (Verner fits are smooth), and the references'
  common OP heritage (Ferguson joins OP near log T 4). The Ferguson-TOPS agreement there is only fair as
  well: Ferguson falls steeply toward its log R = -8 edge at 20 kK while TOPS is flat (`plots/planck_cuts.png`).
- **Rosseland:** within 0.1 dex of both references everywhere in the overlap. It is set by Thomson scattering
  on Saha electrons, plus at most +0.12 dex from lines (largest at 2818 K and 1e-12, where kR/kes = 2.3 with
  Rayleigh and CO).

## 5. Where molecules matter (`plots/molecule_share.png`)
Molecules contribute > 1 % of kP only at log T <= 3.50 and log rho >= -16.2, and > 10 % only at log T 3.45
and log rho >= -14.55. The maximum is 18 % at (3.45, -12), almost entirely CO. Above 3200 K and at
rho < 1e-16 they are irrelevant. H2 and H2O are negligible (< 1e-5 of H, < 1e-11 of O). TiO holds < 3e-6 of
Ti. Caveat: DACE is held at its 2900 K cross sections above 2900 K, an error of a factor of <~ 2 on a <= 18 %
term.

## 6. Comparison with ext2 (`tables/ext2_compare.txt`, `plots/planck_map_vs_ext2.png`)
ext2 is viper's `origin/agcar-opac-1009` (1bf23c1c). The comparison is on the ext2 nodes (log T step 0.025,
log rho step 0.05). Values are this work - ext2 in dex.
- **Planck, log rho -21..-15:** median -0.01, p10/p90 -0.44/+0.50, max 1.21. At rho 1e-20..1e-15 the range is
  -1.02..+0.81.

  | T [K] | rho 1e-20 | 1e-18 | 1e-16 | 1e-15 | (log kP this / ext2) |
  |---|---|---|---|---|---|
  | 2818 | -0.31/-0.34 | -0.28/-0.28 | -0.27/-0.21 | -0.24/-0.18 | |
  | 3981 | 0.73/0.37 | 0.76/0.51 | 0.85/0.66 | 0.93/0.74 | |
  | 5012 | -0.39/0.27 | 1.07/0.61 | 1.51/0.95 | 1.55/1.12 | |
  | 6310 | -0.29/-0.15 | -0.24/0.44 | 0.77/1.03 | 1.50/1.32 | |
  | 7943 | 0.34/-0.47 | 0.57/0.07 | 0.65/0.60 | 0.68/0.87 | |
  | 10000 | 1.14/1.29 | 1.08/1.45 | 1.19/1.61 | 1.31/1.69 | |
  | 12589 | 1.93/1.69 | 1.81/1.87 | 1.81/2.05 | 1.82/2.13 | |
  | 19953 | 1.95/2.96 | 2.51/3.06 | 2.72/3.16 | 2.73/3.21 | |

  The power-law extension (s in [0, 1]) cannot represent the real structure. kP is NOT monotonic in rho at
  fixed T: each time the dominant ion of an abundant element (Mg, Fe, Ca, Si, C) ionises to a closed-shell or
  far-UV ion as rho falls, kP drops by about 1 dex over 2-3 dex in rho. Then it rises or flattens again as the
  next ion takes over (for example 5000 K: 1.07 at 1e-18, -0.39 at 1e-20). ext2 is too high where TOPS was
  clamped (20 kK, up to 1 dex) and too low at 3800-5600 K (0.3-0.65 dex).
- **Rosseland:** median 0.00, |d| <= 0.12 dex. ext2's Saha-electron-scattering extension is confirmed.
- Recommendation for AthenaK: in the low-density corner, replace ext2's Planck values with
  `planck_lowrho_orion_gs98_x0.36_z0.02.txt`, blended at the source floors. Keep ext2's Rosseland as it is.
  Not done here (no table merge was requested).

## 7. Tables (format of TABLES_EXT.md / the AthenaK table reader)
- `tables/planck_lowrho_orion_gs98_x0.36_z0.02.txt` and `tables/rosseland_lowrho_orion_gs98_x0.36_z0.02.txt`.
  Header `# 35 181 3.450 0.025 -21.00 0.05`, then log10 kappa [cm^2/g] with T slowest. The grid is log T
  3.45-4.30 and log rho -21..-12. The nodes coincide with ext2's (`217 421 2.6 0.025 -21 0.05`), so the files
  can be spliced into ext2 without interpolation.
- `tables/planck_components.txt`: kP fractions (lines/bf/ff/H-/molecules), kR, kes, kR without lines, and
  n_e on a coarse subset.
- `tables/sensitivity.txt`, `tables/convergence.txt`, `tables/summary.txt`: all numbers quoted here.
- `scripts/`: `prep_atoms.py` (NIST, gfall, CD-ROM 1, Verner), `merge_lines.py`, `lowrho_lib.py`,
  `chem_check.py` (FastChem), `run_T.py` (one T node), `analyse.py` (tables, validation, ext2, plots),
  `mkcmds.sh` + `run_all.sbatch` (job 213657), and `drv.f` (Verner tabulation driver).
- Rerun recipe:
  1. Download NIST, Kurucz, Verner and Ferguson as listed in the script docstrings.
  2. `prep_atoms.py gf08`, then `prep_atoms.py cd1`, then `merge_lines.py`.
  3. `chem_check.py`.
  4. `sbatch run_all.sbatch`.
  5. `analyse.py`.

## 8. Uncertainty and caveats
- **Numerics:** < 0.001 dex (grid, exact-vs-sampled), plus <= 0.01 dex (wing cut, microturbulence) for kR.
- **Atomic data:** line-list choice <= 0.02 dex. Against the references: about +-0.15 dex at <= 3500 K and
  6300-7100 K, +0.35 dex (Ferguson) at 4000-6000 K, and an unexplained -0.3..-0.55 dex at 11-20 kK. My best
  estimate of the Planck uncertainty is +-0.3 dex below 10 kK, +0.5/-0.1 dex above 11 kK (this table is
  likely low there), and +-1 dex within 9-10 kK, where the H ionisation front makes kP extremely
  T-sensitive.
- **LTE itself:** this is the dominant physical caveat and outside the brief. At rho 1e-20..1e-15 the
  ionisation and excitation are set by the radiation field, not by collisions. The Saha stages used here, for
  example Fe IV and C III at 20 kK and 1e-21, or the 1-dex dips from LTE over-ionisation at 5000-6300 K and
  rho < 1e-18, may not exist in the real atmosphere. The LTE table is a consistent baseline, not a prediction
  of the coupling.
- **Not included:** elements heavier than Zn (Sr, Ba, ... < 0.1 % of lines, by number negligible), grains (none
  at T >= 2800 K), pressure broadening (negligible, section 2), and photoionisation resonances.
