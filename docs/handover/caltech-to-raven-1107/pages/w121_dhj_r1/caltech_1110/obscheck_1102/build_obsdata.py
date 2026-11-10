#!/usr/bin/env python3
# Builds obsdata.json from values transcribed from the papers (see obsdata_notes.md)
# and from data files downloaded into obscheck_1102/lit/ (parsed as plain numbers only).
import json, math, os
import numpy as np

LIT = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/lit'
OUT = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/obsdata.json'

def V(v, e=None, ep=None, em=None):
    d = {'v': v}
    if e is not None: d['e'] = e
    if ep is not None: d['ep'] = ep; d['em'] = em
    return d

def r(x, n=1): return None if x is None else round(float(x), n)

# ---------------------------------------------------------------- stellar params
ST_DELREZ = {'Teff_K': V(6460, 140), 'logg': V(4.2, 0.2), 'FeH': V(0.13, 0.09), 'Rstar_Rsun': V(1.458, 0.030),
             'source': 'Delrez et al. 2016 (MNRAS 458, 4025)'}
ST_SING25 = {'Teff_K': V(6628, 66), 'logg': V(4.251, 0.003), 'FeH': V(0.17, 0.05), 'Rstar_Rsun': V(1.461, 0.005),
             'source': 'Sing et al. 2024 AJ 168, 231 (arXiv 2501.03844) Table 1'}
ST_W76_E20 = {'Teff_K': V(6329, 65), 'logg': V(4.196, 0.106), 'FeH': V(0.366, 0.053), 'Rstar_Rsun': V(1.756, 0.071),
              'source': 'Ehrenreich et al. 2020 Extended Data Table 1'}

W121, W76 = [], []
W121_FE, W76_FE = [], []

# ================================================================ WASP-121b emission
W121.append({
 'id': 'W121_NIRSpec_G395H_MikalEvans2023',
 'citation': {'paper': 'Mikal-Evans et al. 2023, ApJL 943, L17', 'arxiv': '2301.03209', 'doi': '10.3847/2041-8213/acb049',
              'where': 'Sect. 3 text + abstract (white-light, dipole starry map)'},
 'instrument': 'JWST NIRSpec G395H phase curve (GO 1729), first-look analysis', 'status': 'tabulated',
 'bands': [
  {'name': 'NRS1', 'band_um': [2.70, 3.72], 'dayside_ppm': V(3924, 7), 'nightside_ppm': V(136, 8),
   'hotspot_offset_deg_east': V(3.36, 0.11), 'Tday_K': V(2762, ep=30, em=32), 'Tnight_K': V(926, 12), 'RpRs': V(0.122551, 0.000063)},
  {'name': 'NRS2', 'band_um': [3.82, 5.15], 'dayside_ppm': V(4924, 9), 'nightside_ppm': V(630, 10),
   'hotspot_offset_deg_east': V(2.66, 0.12), 'Tday_K': V(2768, 39), 'Tnight_K': V(1122, 10), 'RpRs': V(0.123668, 0.000068)}],
 'offset_convention_original': 'phase-curve peak this many degrees BEFORE mid-eclipse (eastward); same as ours',
 'stellar_params': dict(ST_DELREZ, logg=V(4.2, 0.2), note='Tb uses ATLAS9 stellar models drawn from Delrez 2016 posteriors'),
 'notes': 'Superseded by Evans-Soma 2025 full analysis (offsets 2.95/2.39 deg); peak-to-trough = dayside - nightside (dipole map).'})

# --- Evans-Soma 2025: Zenodo phase-resolved spectra
Z = os.path.join(LIT, 'zenodo_15377434')
wav = np.loadtxt(os.path.join(Z, 'emspec_wav.txt'))            # (349,2) edges
phi = np.loadtxt(os.path.join(Z, 'emspec_phi.txt'))            # (36,3) id, lo, hi ; 0=transit, 0.5=eclipse
fp = np.loadtxt(os.path.join(Z, 'emspec_fpfs_vals.txt'))       # (36,349)
fe = np.loadtxt(os.path.join(Z, 'emspec_fpfs_uncs.txt'))
tb = np.loadtxt(os.path.join(Z, 'emspec_tbright_vals.txt'))
tbl = np.loadtxt(os.path.join(Z, 'emspec_tbright_uncs_lower.txt'))
tbu = np.loadtxt(os.path.join(Z, 'emspec_tbright_uncs_upper.txt'))
assert fp.shape == (36, wav.shape[0])
wc = wav.mean(1); hw = 0.5 * (wav[:, 1] - wav[:, 0])
nrs1 = wc < 3.76; nrs2 = wc > 3.78
ph_lo = np.mod(phi[:, 1], 1.0); ph_hi = np.mod(phi[:, 2], 1.0); ph_hi[ph_hi == 0] = 1.0
ph_c = np.mod(0.5 * (phi[:, 1] + phi[:, 2]), 1.0)
def bandavg(row, erow, m):
    n = m.sum(); return float(row[m].mean()), float(math.sqrt((erow[m] ** 2).sum()) / n)
pc = []
for k in range(36):
    b1 = bandavg(fp[k], fe[k], nrs1); b2 = bandavg(fp[k], fe[k], nrs2)
    pc.append({'n_phi': int(phi[k, 0]), 'phase_center': r(ph_c[k], 4), 'phase_lo': r(ph_lo[k], 4), 'phase_hi': r(ph_hi[k], 4),
               'NRS1_ppm': V(r(b1[0]), r(b1[1])), 'NRS2_ppm': V(r(b2[0]), r(b2[1]))})
def spec(k):
    return {'n_phi': int(phi[k, 0]), 'phase_range': [r(ph_lo[k], 4), r(ph_hi[k], 4)],
            'wl_um': [r(x, 4) for x in wc], 'halfwidth_um': [r(x, 5) for x in hw],
            'fpfs_ppm': [r(x) for x in fp[k]], 'err_ppm': [r(x) for x in fe[k]],
            'Tb_K': [r(x, 0) for x in tb[k]], 'Tb_err_lo_K': [r(x, 0) for x in tbl[k]], 'Tb_err_hi_K': [r(x, 0) for x in tbu[k]]}
kday = int(np.where(phi[:, 0] == 34)[0][0]); knight = int(np.where(phi[:, 0] == 18)[0][0])
W121.append({
 'id': 'W121_NIRSpec_G395H_EvansSoma2025',
 'citation': {'paper': 'Evans-Soma et al. 2025, Nature Astronomy (SiO and super-stellar C/O in WASP-121b)', 'arxiv': '2506.01771',
              'where': 'Fig. 1 (n_phi=34 dayside, n_phi=18 nightside), Extended Data Fig. 2; data = Zenodo 10.5281/zenodo.15377434 (emspec_*.txt)'},
 'instrument': 'JWST NIRSpec G395H phase curve (GO 1729), full reanalysis', 'status': 'tabulated (Zenodo); band averages derived here',
 'white_light_offsets_deg_east': {'NRS1': V(2.95, ep=0.11, em=0.12), 'NRS2': V(2.39, ep=0.12, em=0.13),
     'source': 'Supplementary Table 1 of Evans-Soma 2025 as quoted by Splinter et al. 2025 (arXiv 2509.09760) Sect. 5; spectroscopic offsets all < 10 deg'},
 'dayside_spectrum': spec(kday), 'nightside_spectrum': spec(knight),
 'band_avg_phase_curve': {'definition': 'unweighted mean over channels of each detector (2.725-3.72 / 3.83-5.17 um), err = sqrt(sum err^2)/N; derived here from Zenodo arrays',
                          'points': pc},
 'phase_resolved_spectra_file': os.path.join(Z, 'emspec_fpfs_vals.txt') + ' (36 phase bins x 349 channels; rows=phase bins in emspec_phi.txt order)',
 'stellar_params': {'source': 'Supplementary Table 9 of the paper (not in arXiv text); Kahle 2026 states they are compatible with Delrez 2016 (Teff 6459+-140)'},
 'notes': 'Phase bins 46 min wide; bins adjacent to eclipse (n_phi 34/35) include the partial-ingress-free dayside. Phase in Zenodo file: 0 transit, +-0.5 eclipse; converted to [0,1).'})

# --- Splinter 2025 NIRISS SOSS
dig_day = [(0.65, 300), (0.75, 400), (0.85, 550), (0.95, 680), (1.05, 760), (1.15, 900), (1.25, 1000), (1.35, 1150), (1.45, 1350),
           (1.55, 1380), (1.65, 1450), (1.75, 1620), (1.85, 2000), (1.95, 2050), (2.05, 2120), (2.15, 2230), (2.25, 2200),
           (2.35, 2650), (2.45, 2950), (2.55, 3150), (2.65, 3250), (2.75, 3300)]
dig_night_exo = [(0.7, 20), (0.9, 0), (1.1, 20), (1.3, 0), (1.5, 20), (1.7, 60), (1.9, 100), (2.1, 120), (2.3, 150), (2.5, 350), (2.7, 400)]
dig_night_nam = [(0.7, -250), (0.9, -200), (1.1, -150), (1.3, -120), (1.5, -100), (1.7, -80), (1.9, -20), (2.1, 0), (2.3, 0), (2.5, 50), (2.7, 150)]
W121.append({
 'id': 'W121_NIRISS_SOSS_Splinter2025',
 'citation': {'paper': 'Splinter et al. 2025, AJ 170, 323', 'arxiv': '2509.09760', 'doi': '10.3847/1538-3881/ae0e52',
              'where': 'Table 1 (white light), Sect. 5 text (offsets, Teff), Fig. Nightside-v-Dayside (spectra, digitised)'},
 'instrument': 'JWST NIRISS SOSS phase curve (GTO 1201)', 'status': 'tabulated (white light); spectra digitised',
 'bands': [
  {'name': 'Order 1 white', 'band_um': [0.85, 2.85], 'dayside_ppm': V(1155.99, ep=14.30, em=14.04),
   'hotspot_offset_deg_east': V(5.1, 1.4), 'RpRs': V(0.1222, 0.00014)},
  {'name': 'Order 2 white', 'band_um': [0.60, 0.85], 'dayside_ppm': V(363.13, ep=20.79, em=21.64),
   'hotspot_offset_deg_east': V(10.5, 9.9), 'RpRs': V(0.12293, 0.0002)}],
 'bolometric': {'Tday_eff_K': V(2717, 17), 'Tnight_eff_K': V(1562, ep=18, em=19), 'Bond_albedo': V(0.277, 0.016), 'epsilon': V(0.246, 0.014),
                'geometric_albedo_0.6-0.85um': V(0.093, ep=0.029, em=0.027), 'EBM_best': 'A_B~0.31, wind ~0.2 km/s, P0~1 bar'},
 'offset_convention_original': 'positive = eastward (max before mid-eclipse); same as ours',
 'dayside_spectrum_digitised': {'wl_um': [a for a, b in dig_day], 'halfwidth_um': 0.05, 'fpfs_ppm': [b for a, b in dig_day], 'err_ppm': 120,
                                'note': 'eye-averaged 0.1-um means of the exoTEDRF+NAMELESS points; digitisation error ~100-150 ppm. Better: use Kahle 2026 Zenodo NIRISS eclipse spectrum (entry W121_JWST_panchromatic_Kahle2026).'},
 'nightside_spectrum_digitised': {'exoTEDRF': {'wl_um': [a for a, b in dig_night_exo], 'fpfs_ppm': [b for a, b in dig_night_exo]},
                                  'NAMELESS': {'wl_um': [a for a, b in dig_night_nam], 'fpfs_ppm': [b for a, b in dig_night_nam]},
                                  'err_ppm': 100, 'note': 'reduction-dependent (exoTEDRF penalises negative flux; NAMELESS allows it); treat as ~0-400 ppm, rising beyond 2.3 um'},
 'stellar_params': dict(ST_DELREZ, note='Teff 6459+-140 (Delrez) used for Tb; a/R*=3.8002, i=88.345'),
 'notes': 'Nightside white-light flux not tabulated. Offsets near zero for lambda>1.5 um; larger eastward offsets at short lambda.'})

# --- Kahle 2026 panchromatic dayside (Zenodo 20767846)
K = os.path.join(LIT, 'zenodo_20767846', 'WASP-121b_dayside_spectra')
def rdcsv(fn):
    rows = []
    for line in open(os.path.join(K, fn)):
        if line.startswith('#') or not line.strip(): continue
        rows.append([float(x) for x in line.strip().split(',')])
    return np.array(rows)
miri = rdcsv('WASP121b_emission_MIRI_LRS.txt'); nsp = rdcsv('WASP121b_emission_NIRSpec_G395H.txt')
n1 = rdcsv('WASP121b_emission_niriss_firstorder.txt'); n2 = rdcsv('WASP121b_emission_niriss_secondorder.txt')
def mk(a, hw=None):
    if a.shape[1] == 4:   # wl, binwidth, fpfs, err
        return {'wl_um': [r(x, 4) for x in a[:, 0]], 'halfwidth_um': [r(0.5 * x, 5) for x in a[:, 1]],
                'fpfs_ppm': [r(1e6 * x) for x in a[:, 2]], 'err_ppm': [r(1e6 * x) for x in a[:, 3]]}
    d = {'wl_um': [r(x, 4) for x in a[:, 0]], 'fpfs_ppm': [r(1e6 * x) for x in a[:, 1]], 'err_ppm': [r(1e6 * x) for x in a[:, 2]]}
    dw = np.gradient(a[:, 0]); d['halfwidth_um'] = [r(0.5 * x, 5) for x in dw]
    return d
W121.append({
 'id': 'W121_JWST_panchromatic_Kahle2026',
 'citation': {'paper': 'Kahle et al. 2026 (arXiv preprint, Jun 2026) "The panchromatic JWST dayside spectrum of WASP-121 b"', 'arxiv': '2606.21855',
              'where': 'Table 1 (MIRI white light), Fig. 2; data = Zenodo 10.5281/zenodo.20767846 (WASP-121b_dayside_spectra.zip)'},
 'instrument': 'JWST MIRI LRS eclipse (GO 2961) + NIRISS SOSS eclipse (Pelletier 2026 NAMELESS reduction) + NIRSpec G395H (Evans-Soma 2025 phase bins converted to eclipse phase)',
 'status': 'tabulated (Zenodo)',
 'bands': [{'name': 'MIRI LRS white', 'band_um': [5.0, 12.0], 'dayside_ppm': V(6188, 36)}],
 'dayside_spectrum': {'MIRI_LRS': mk(miri), 'NIRSpec_G395H': mk(nsp), 'NIRISS_order1': mk(n1), 'NIRISS_order2': mk(n2)},
 'stellar_params': {'Teff_K': V(6459, 140), 'logg': V(4.251, 0.003), 'FeH': V(0.17, 0.05), 'Rstar_Rsun': V(1.461, 0.005),
                    'source': 'Kahle 2026 Table 2 (Delrez Teff; Sing 2024 logg, [M/H], R*)', 'RpRs_MIRI_fixed': 0.123180},
 'notes': 'NIRISS/MIRI halfwidths from channel spacing (derived); NIRSpec column 2 = bin width. Retrieval found NIRSpec offset 260+-120 ppm and MIRI 237+-75 ppm vs NIRISS anchor (not preferred by model selection).'})

# --- Sing 2025 planet RV from NIRSpec phase curve (dayside emission lines)
W121.append({
 'id': 'W121_NIRSpec_RV_Sing2025',
 'citation': {'paper': 'Sing et al. 2024, AJ 168, 231', 'arxiv': '2501.03844', 'doi': '10.3847/1538-3881/ad7fe7', 'where': 'Table 1, Sect. 5'},
 'instrument': 'JWST NIRSpec G395H phase curve, planet emission-line RV', 'status': 'tabulated',
 'Kp_kms': {'combined': V(215.7, 1.1), 'NRS1': V(212.7, 1.8), 'NRS2': V(218.2, 1.0), 'expected_orbital': V(217.8, 1.0)},
 'notes': 'NRS1 lower than NRS2 by 5.5+-2.2 km/s; authors attribute ~1 km/s to rotation, 4.5+-2.2 to winds (de-projected ~5 km/s zonal). Also gives stellar params (Teff 6628+-66, logg 4.251, [Fe/H] 0.17, R* 1.461, Mp 1.170+-0.043 MJ, Rp 1.742 RJ).',
 'stellar_params': ST_SING25})

# --- HST WFC3 Mikal-Evans 2022 (Nature Astronomy source data)
M = os.path.join(LIT, 'nature_me22')
def rd(fn): return np.loadtxt(os.path.join(M, fn))
dfix = rd('m4/planetEmission_c2Fixed_Dayside.txt'); nfix = rd('m4/planetEmission_c2Fixed_Nightside.txt')
dfree = rd('m4/planetEmission_c2Free_Dayside.txt')
bb = rd('m3/broadbandLightCurve_jointBinned.txt')
def sp4(a): return {'wl_um': [r(x, 4) for x in a[:, 0]], 'halfwidth_um': [r(0.5 * x, 5) for x in a[:, 1]],
                    'fpfs_ppm': [r(x) for x in a[:, 2]], 'err_ppm': [r(x) for x in a[:, 3]]}
phase_specs = []
for fn in sorted(os.listdir(os.path.join(M, 'm7'))):
    if not fn.startswith('planetEmission_Binned_phase'): continue
    a = rd('m7/' + fn); p = float(fn.split('phase')[1].replace('.txt', ''))
    d = sp4(a); d['phase'] = p; phase_specs.append(d)
tb_me22 = {0.05: (1555, 52), 0.12: (1972, 22), 0.17: (2190, 16), 0.23: (2366, 13), 0.28: (2505, 11), 0.32: (2622, 10), 0.38: (2699, 9),
           0.43: (2749, 8), 0.57: (2670, 9), 0.62: (2580, 9), 0.68: (2461, 12), 0.72: (2331, 13), 0.78: (2178, 16), 0.82: (2022, 20),
           0.88: (1838, 28), 0.95: (1651, 38)}
bbpts = []
for p, v, e in bb:
    pm = p % 1.0
    ev = (min(pm, 1 - pm) < 0.06) or (abs(pm - 0.5) < 0.06)
    bbpts.append({'phase': r(pm, 4), 'fp_ppm': r(v), 'err_ppm': r(e), 'in_or_near_transit_eclipse': bool(ev)})
W121.append({
 'id': 'W121_HST_WFC3_G141_MikalEvans2022',
 'citation': {'paper': 'Mikal-Evans et al. 2022, Nature Astronomy 6, 471', 'arxiv': '2202.09884', 'doi': '10.1038/s41550-021-01592-w',
              'where': 'Extended Data Table 1 (broadband fit), Source Data zips (MOESM3/4/7), Supplementary Table 6 (Tb vs phase)'},
 'instrument': 'HST WFC3 G141, two full phase curves (2018, 2019)', 'status': 'tabulated (Source Data)',
 'bands': [{'name': 'G141 broadband', 'band_um': [1.12, 1.64],
            'hotspot_offset_deg_east': V(6.4, ep=1.7, em=2.3),
            'fit_params_2nd_order': {'c1_ppm': V(1158, ep=26, em=30), 'c2_rad': V(3.03, ep=0.03, em=0.04), 'RpRs': 0.12199, 'AF_percent': V(95.1, 2.5)},
            'note': 'offset = pi - c2 (c2 = phase of peak, pi at eclipse); preferred 2nd-order sinusoid. 1st-order fit: c2=3.06 rad -> 4.7 deg east. Errors on c1/c2 partly garbled in arXiv PDF; offset error from Splinter 2025 quote.'}],
 'dayside_spectrum': dict(sp4(dfix), note='c2 fixed; dayside = average of emission just before/after eclipse'),
 'dayside_spectrum_c2free': dict(sp4(dfree), note='includes TESS-band point at 0.7865 um (409+-21 ppm)'),
 'nightside_spectrum': dict(sp4(nfix), phase=0.95, note='phase 0.90-1.00 bin, dayside crescent < 8 ppm'),
 'phase_resolved_spectra': phase_specs,
 'Tb_vs_phase_K': [{'phase': k, 'Tb': V(v, e)} for k, (v, e) in sorted(tb_me22.items())],
 'broadband_phase_curve_binned': {'points': bbpts, 'note': 'planet emission (ppm) from joint 2018+2019 binned light curve; values in/near transit (negative) and eclipse are not planet emission'},
 'stellar_params': {'source': 'not extracted (see Methods)'},
 'notes': 'Changeat et al. 2024 ApJS 270, 34 (arXiv 2401.01465) reanalysed the same data + 2016 transit/eclipses with a 3-region model: best evidence for fixed hotspot offset 30 deg east, size 50 deg; no tabulated spectra; reports epoch-to-epoch differences.'})

W121.append({
 'id': 'W121_HST_WFC3_G141_eclipse_MikalEvans2020',
 'citation': {'paper': 'Mikal-Evans et al. 2020, MNRAS 496, 1638', 'arxiv': '2005.09631', 'doi': '10.1093/mnras/staa1628', 'where': 'Table 2 (spectroscopic), Table 1 (white)'},
 'instrument': 'HST WFC3 G141, 5 eclipses (incl. Evans 2017 visit)', 'status': 'tabulated',
 'bands': [{'name': 'G141 white', 'band_um': [1.12, 1.64], 'dayside_ppm': V(1150, ep=21, em=19)}],
 'dayside_spectrum': (lambda rows: {'wl_um': [r(0.5 * (a + b), 4) for a, b, c, d, e in rows], 'halfwidth_um': [r(0.5 * (b - a), 4) for a, b, c, d, e in rows],
                                    'fpfs_ppm': [c for a, b, c, d, e in rows], 'err_ppm_plus': [d for a, b, c, d, e in rows], 'err_ppm_minus': [e for a, b, c, d, e in rows]})(
   [(1.120, 1.138, 903, 53, 52), (1.138, 1.157, 991, 59, 60), (1.157, 1.175, 1002, 58, 56), (1.175, 1.194, 1029, 50, 49), (1.194, 1.212, 1066, 58, 58),
    (1.212, 1.231, 983, 54, 57), (1.231, 1.249, 1031, 48, 51), (1.249, 1.268, 1015, 55, 50), (1.268, 1.286, 994, 49, 48), (1.286, 1.305, 1028, 55, 53),
    (1.305, 1.323, 1008, 55, 60), (1.323, 1.342, 1077, 55, 56), (1.342, 1.360, 1160, 51, 52), (1.360, 1.379, 1110, 61, 57), (1.379, 1.397, 1262, 57, 59),
    (1.397, 1.416, 1360, 57, 56), (1.416, 1.434, 1193, 54, 56), (1.434, 1.453, 1304, 51, 54), (1.453, 1.471, 1331, 58, 63), (1.471, 1.490, 1342, 57, 57),
    (1.490, 1.508, 1304, 62, 60), (1.508, 1.527, 1276, 62, 62), (1.527, 1.545, 1210, 66, 65), (1.545, 1.564, 1307, 62, 60), (1.564, 1.582, 1388, 63, 66),
    (1.582, 1.601, 1299, 69, 69), (1.601, 1.619, 1270, 64, 63), (1.619, 1.638, 1286, 68, 68)]),
 'notes': 'Supersedes Evans et al. 2017 Nature 548, 58 (arXiv 1708.01076) single-eclipse G141 spectrum (not transcribed).'})

# --- Spitzer
W121.append({
 'id': 'W121_Spitzer_Morello2023',
 'citation': {'paper': 'Morello et al. 2023, A&A 676, A54', 'arxiv': '2307.00669', 'doi': '10.1051/0004-6361/202346643', 'where': 'Table 4 (posteriors), Table 1 (system)'},
 'instrument': 'Spitzer IRAC 3.6/4.5 um phase curves (PID 13242)', 'status': 'tabulated',
 'bands': [
  {'name': 'IRAC1', 'band_um': [3.16, 3.92], 'dayside_max_ppm': V(4230, 80), 'nightside_min_ppm': V(50, 240),
   'F_phase0.5_ppm': V(4210, 80), 'F_phase0_ppm': V(220, 260), 'hotspot_offset_deg_east': V(-5.9, 1.6), 'min_offset_deg_before_transit': V(31, 6),
   'Tday_K': V(2670, ep=55, em=40), 'Tnight_K': V(710, ep=270, em=710), 'Tb_phase0_K': V(810, ep=360, em=810), 'RpRs': V(0.1228, 0.0005)},
  {'name': 'IRAC2', 'band_um': [3.95, 5.06], 'dayside_max_ppm': V(5090, 90), 'nightside_min_ppm': V(710, 270),
   'F_phase0.5_ppm': V(5080, 90), 'F_phase0_ppm': V(720, 270), 'hotspot_offset_deg_east': V(-5.0, ep=3.1, em=3.4), 'min_offset_deg_before_transit': V(-4.5, 2.9),
   'Tday_K': V(2700, ep=70, em=50), 'Tnight_K': V(1130, ep=130, em=160), 'RpRs': V(0.1231, 0.0005)}],
 'offset_convention_original': 'Delta phi positive = AFTER mid-eclipse (westward peak) / after mid-transit for the minimum; sign flipped here',
 'stellar_params': ST_DELREZ, 'notes': 'Band edges = nominal IRAC half-power ranges (not from paper). Westward offsets disputed by Davenport 2025 and JWST.'})
W121.append({
 'id': 'W121_Spitzer_Davenport2025',
 'citation': {'paper': 'Davenport et al. 2025 (arXiv preprint)', 'arxiv': '2503.12521', 'where': 'Table 3 (Fit Results), Sect. 4.1 text'},
 'instrument': 'Spitzer IRAC 3.6/4.5 um phase curves (same data as Morello 2023)', 'status': 'tabulated',
 'bands': [
  {'name': 'IRAC1', 'band_um': [3.16, 3.92], 'dayside_ppm': V(4074, 60), 'eclipse_depth_ppm': V(4077, 59), 'nightside_ppm': V(554, 96),
   'max_ppm': V(4085, 96), 'min_ppm': V(553, 96), 'semi_amplitude_ppm': V(1771, 95), 'hotspot_offset_deg_east': V(-0.78, 1.87),
   'Tday_K_table': V(2779, 40), 'Tnight_K_table': V(1259, 67), 'Tday_K_text': V(2543, 42), 'Tnight_K_text': V(1082, 61), 'RpRs': V(0.12284, 0.00081)},
  {'name': 'IRAC2', 'band_um': [3.95, 5.06], 'dayside_ppm': V(5116, 75), 'eclipse_depth_ppm': V(5121, 76), 'nightside_ppm': V(1046, 107),
   'max_ppm': V(5141, 107), 'min_ppm': V(1045, 107), 'semi_amplitude_ppm': V(2048, 109), 'hotspot_offset_deg_east': V(0.42, 1.74),
   'Tday_K_table': V(2905, 51), 'Tnight_K_table': V(1349, 54), 'RpRs': V(0.12402, 0.00042)}],
 'offset_convention_original': 'positive = eastward',
 'stellar_params': {'Teff_K': V(6776, 138), 'logg': V(4.234, 0.11), 'FeH': V(0.13, 0.09), 'Rstar_Rsun': V(0.910, ep=0.025, em=0.023),
                    'source': 'Davenport Table 1 AS PRINTED (R*=0.91 looks like a typo/swap; Teff unusual)'},
 'notes': 'Amplitude = (max-min)/2. Tday/Tnight in table disagree with text (2543/1082 K at 3.6 um) - use fluxes, not their T.'})
W121.append({
 'id': 'W121_Spitzer45_Dang2025',
 'citation': {'paper': 'Dang et al. 2025, AJ (Comprehensive Analysis of Spitzer 4.5 um Phase Curves)', 'arxiv': '2408.13308', 'where': 'Table 3 (best fit), Table 4 (derived), Table 6 (system)'},
 'instrument': 'Spitzer IRAC 4.5 um (BLISS detrending)', 'status': 'tabulated',
 'bands': [{'name': 'IRAC2', 'band_um': [3.95, 5.06], 'eclipse_depth_ppm': V(4800, 100), 'min_ppm': V(600, ep=400, em=300), 'max_ppm': V(4900, 100),
            'hotspot_offset_deg_east': V(-5, 4), 'Tday_K': V(2680, 60), 'Tnight_K': V(1100, 200), 'Bond_albedo': V(0.31, ep=0.06, em=0.07),
            'epsilon': V(0.07, ep=0.07, em=0.04), 'RpRs': V(0.1249, 0.0009)}],
 'offset_convention_original': 'degrees East (positive = east)',
 'stellar_params': {'Teff_K': V(6500, 100), 'Rstar_Rsun': V(1.46, 0.03), 'source': 'Dang 2025 Table 6'}})

# --- TESS
W121.append({
 'id': 'W121_TESS_Daylan2021',
 'citation': {'paper': 'Daylan et al. 2021, AJ 161, 131', 'arxiv': '1909.03000', 'doi': '10.3847/1538-3881/abd8d2', 'where': 'Table 3 (nominal model column) + macros'},
 'instrument': 'TESS phase curve (sector 7)', 'status': 'tabulated',
 'bands': [{'name': 'TESS', 'band_um': [0.6, 0.95], 'eclipse_depth_ppm': V(482, ep=39, em=41), 'nightside_ppm': V(65, ep=40, em=37),
            'atm_amplitude_ppm': V(418, 34), 'hotspot_offset_deg_east': V(-6.3, 3.9), 'Tday_K': V(3012, ep=40, em=42),
            'Tnight_K': V(2022, ep=254, em=602), 'geometric_albedo': V(0.07, ep=0.037, em=0.040)}],
 'offset_convention_original': 'reported as westward 6.3+-3.9 deg; sign flipped', 'stellar_params': {'source': 'Delrez 2016 priors'}})
W121.append({
 'id': 'W121_TESS_Bourrier2020',
 'citation': {'paper': 'Bourrier et al. 2020, A&A 637, A36', 'arxiv': '1909.03010', 'doi': '10.1051/0004-6361/201936647', 'where': 'Table 1'},
 'instrument': 'TESS phase curve', 'status': 'tabulated',
 'bands': [{'name': 'TESS', 'band_um': [0.6, 1.0], 'dayside_ppm': V(419, ep=47, em=41), 'nightside_ppm': V(63, 27),
            'hotspot_offset_deg_east': 'consistent with 0 (xi = -0.016+0.061-0.064)', 'Tday_K': V(2870, 50), 'Tnight_K': V(2190, ep=106, em=145),
            'RpRs': V(0.12355, ep=0.00033, em=0.00029)}],
 'stellar_params': {'source': 'Delrez 2016 (M*, R*)'}})

# ================================================================ WASP-76b emission
W76.append({
 'id': 'W76_Spitzer_MayStevenson2021',
 'citation': {'paper': 'May, Komacek, Stevenson et al. 2021, AJ 162, 158', 'arxiv': '2107.03349', 'doi': '10.3847/1538-3881/ac0e30',
              'where': 'Sect. 3.2-3.3 text, Table 3 (eclipses), Table A1 (fit params)'},
 'instrument': 'Spitzer IRAC 3.6 (2 phase curves, joint) / 4.5 um (1 phase curve)', 'status': 'tabulated; nightside flux derived',
 'bands': [
  {'name': 'IRAC1', 'band_um': [3.16, 3.92], 'eclipse_depth_ppm': V(2539, 30), 'phase_amplitude_ppm': V(804.0, 42.5),
   'nightside_ppm_derived': V(931, 96), 'hotspot_offset_deg_east': V(0.68, 0.48), 'Tday_K': V(2471, 27), 'Tnight_K': V(1518, 61)},
  {'name': 'IRAC2', 'band_um': [3.95, 5.06], 'eclipse_depth_ppm': V(3729, 52), 'phase_amplitude_ppm': V(1464, 38),
   'nightside_ppm_derived': V(801, 89), 'hotspot_offset_deg_east': V(0.67, 0.2), 'Tday_K': V(2699, 32), 'Tnight_K': V(1259, 44)}],
 'offset_convention_original': 'SPIDERMAN hotspot longitude lo0 (positive assumed east; paper calls offsets small, no sign discussion)',
 'stellar_params': {'Teff_K': 6250.0, 'source': 'West et al. 2016 (Kurucz model for Tb; Ts fixed 6250 in fit table)', 'RpRs_fit': [0.106, 0.10048]},
 'notes': 'Dilution (companion) corrected. Nightside derived as eclipse - 2*amplitude assuming amplitude = semi-amplitude (gives ~930/800 ppm, consistent with Dang 2025 Fmin 900+-300); verify before strict use. Authors flag the 3.6 um data as less reliable.'})
W76.append({
 'id': 'W76_Spitzer45_Dang2025',
 'citation': {'paper': 'Dang et al. 2025, AJ', 'arxiv': '2408.13308', 'where': 'Tables 3, 4, 6'},
 'instrument': 'Spitzer IRAC 4.5 um (Poly5 detrending)', 'status': 'tabulated',
 'bands': [{'name': 'IRAC2', 'band_um': [3.95, 5.06], 'eclipse_depth_ppm': V(3400, 100), 'min_ppm': V(900, 300), 'max_ppm': V(3400, 100),
            'hotspot_offset_deg_east': V(-3, ep=7, em=9), 'Tday_K': V(2580, 60), 'Tnight_K': V(1400, 200), 'Bond_albedo': V(0.16, ep=0.09, em=0.1),
            'epsilon': V(0.2, 0.1), 'RpRs': V(0.1049, ep=0.0008, em=0.0009)}],
 'stellar_params': {'Teff_K': V(6200, 100), 'Rstar_Rsun': V(1.73, 0.04), 'source': 'Dang 2025 Table 6'},
 'notes': 'Eclipse depth likely NOT dilution-corrected for WASP-76B (3400 vs 3729 in May 2021) - check.'})
fu_rows = [(1.1518, 490, 49), (1.1704, 573, 49), (1.1890, 570, 48), (1.2076, 540, 48), (1.2262, 553, 47), (1.2448, 550, 47), (1.2634, 458, 53), (1.2820, 477, 50),
           (1.3006, 567, 47), (1.3192, 659, 45), (1.3378, 640, 47), (1.3564, 665, 47), (1.3750, 740, 46), (1.3936, 716, 48), (1.4122, 746, 47), (1.4308, 699, 49),
           (1.4494, 789, 47), (1.4680, 820, 57), (1.4866, 868, 48), (1.5052, 947, 51), (1.5238, 874, 50), (1.5424, 912, 59), (1.5610, 830, 52), (1.5796, 906, 52),
           (1.5982, 911, 57), (1.6168, 943, 56)]
W76.append({
 'id': 'W76_HST_WFC3_Fu2021',
 'citation': {'paper': 'Fu et al. 2021, AJ 162, 108 (Hubble PanCET WASP-76b)', 'arxiv': '2005.02568', 'doi': '10.3847/1538-3881/ac1200', 'where': 'Table 7 (emission spectrum), Table 3 (Spitzer)'},
 'instrument': 'HST WFC3 G141 eclipse + Spitzer IRAC eclipses', 'status': 'tabulated',
 'dayside_spectrum': {'wl_um': [a for a, b, c in fu_rows], 'halfwidth_um': 0.0093, 'fpfs_ppm': [b for a, b, c in fu_rows], 'err_ppm': [c for a, b, c in fu_rows],
                      'note': 'dilution-corrected (table also lists dilution factors 1.048-1.066)'},
 'bands': [{'name': 'IRAC1 (3 eclipses, diluted-corrected)', 'band_um': [3.16, 3.92], 'eclipse_depth_ppm_each': [V(2827, 69), V(3082, 102), V(3299, 94)]},
           {'name': 'IRAC2', 'band_um': [3.95, 5.06], 'eclipse_depth_ppm': V(3665, 89)}],
 'stellar_params': {'Teff_K': V(6366, ep=92, em=90), 'logg': V(4.122, ep=0.016, em=0.020), 'FeH': V(0.204, ep=0.093, em=0.096), 'Rstar_Rsun': V(1.744, ep=0.045, em=0.042),
                    'RpRs': V(0.10873, ep=0.00048, em=0.00047), 'source': 'Fu 2021 EXOFASTv2 table'}})
ed_rows = [(1.12625, 0.0219, 607, 39), (1.14775, 0.0211, 711, 40), (1.16860, 0.0206, 736, 43), (1.18880, 0.0198, 597, 48), (1.20835, 0.0193, 772, 53),
           (1.22750, 0.0190, 544, 41), (1.24645, 0.0189, 633, 50), (1.26550, 0.0192, 692, 42), (1.28475, 0.0193, 742, 47), (1.30380, 0.0188, 770, 53),
           (1.32260, 0.0188, 836, 52), (1.34145, 0.0189, 870, 50), (1.36050, 0.0192, 976, 52), (1.38005, 0.0199, 903, 44), (1.40000, 0.0200, 882, 48),
           (1.42015, 0.0203, 955, 51), (1.44060, 0.0206, 988, 48), (1.46150, 0.0212, 1139, 45), (1.48310, 0.0220, 1160, 44), (1.50530, 0.0224, 1035, 51),
           (1.52800, 0.0230, 1067, 45), (1.55155, 0.0241, 1038, 48), (1.57625, 0.0253, 1138, 53), (1.60210, 0.0264, 1067, 56), (1.62945, 0.0283, 1031, 54)]
W76.append({
 'id': 'W76_HST_WFC3_Edwards2020',
 'citation': {'paper': 'Edwards et al. 2020, AJ 160, 8 (ARES I)', 'arxiv': '2005.02374', 'doi': '10.3847/1538-3881/ab9225', 'where': 'Table 3 (corrected spectra)'},
 'instrument': 'HST WFC3 G141 eclipse (same visit as Fu 2021)', 'status': 'tabulated',
 'dayside_spectrum': {'wl_um': [a for a, b, c, d in ed_rows], 'halfwidth_um': [round(b / 2, 5) for a, b, c, d in ed_rows],
                      'fpfs_ppm': [c for a, b, c, d in ed_rows], 'err_ppm': [d for a, b, c, d in ed_rows],
                      'note': 'Depth[%] x 1e4; corrected for companion dilution with their (larger, 1.08-1.11) correction factor -> ~10-20% higher than Fu 2021'},
 'stellar_params': ST_W76_E20})
W76.append({'id': 'W76_JWST_none', 'status': 'not found',
            'notes': 'No published JWST eclipse/phase curve of WASP-76b found as of 2026-11-02 (NIRSpec G395H phase curve = Cycle 3 GO 5268, MIRI LRS transit GO 6045, results unpublished). CHEOPS/TESS occultations exist (A&A 684, A27, 2024: TESS 260+-11 ppm, CHEOPS 152+-10 ppm; not checked in detail).'})

# ================================================================ Fe I transmission, WASP-76b
P76_h = 1.80988198 * 24.0
e20_v = [(-1.71, -1.2), (-1.52, -4.9), (-1.31, -2.6), (-1.10, -5.0), (-0.89, -5.8), (-0.69, -7.1), (-0.46, -9.4), (-0.26, -9.7),
         (0.67, -9.4), (0.88, -9.4), (1.09, -10.6), (1.29, -11.0), (1.50, -11.3), (1.70, -11.5)]
e20_d = [(-1.71, 300), (-1.52, 600), (-1.31, 340), (-1.10, 500), (-0.89, 600), (-0.69, 500), (-0.46, 400), (-0.26, 507),
         (0.67, 590), (0.88, 605), (1.09, 600), (1.29, 775), (1.50, 705), (1.70, 565)]
W76_FE.append({
 'id': 'W76_FeCCF_Ehrenreich2020',
 'citation': {'paper': 'Ehrenreich et al. 2020, Nature 580, 597', 'arxiv': '2003.05528', 'doi': '10.1038/s41586-020-2107-1', 'where': 'Extended Data Fig. 7 (digitised black binned curve), main text'},
 'instrument': 'ESPRESSO 1-UT, 2 transits; stellar-mask CCF (F9 mask, mostly Fe I)', 'Kp_kms': V(196.52, 0.94), 'vsys_kms': -1.11,
 'status': 'tabulated summary + digitised time series',
 'summary': {'ingress_velocity_kms': '0 to -5', 'plateau_velocity_kms': V(-11, 0.7), 'plateau_from_h': -0.4,
             'mean_depth_ppm': V(494, 27), 'first_half_depth_ppm': V(434, 32), 'end_depth_ppm': V(628, 49), 'FWHM_kms': V(8.6, 0.7)},
 'time_series_digitised': {'t_from_midtransit_h': [a for a, b in e20_v], 'phase': [round(a / P76_h, 4) for a, b in e20_v],
                           'v_kms': [b for a, b in e20_v], 'v_err_kms': 1.0,
                           'depth_ppm': [b for a, b in e20_d], 'depth_err_ppm': 100,
                           'note': 'black curve = both epochs binned by 2; grey band ~+-0.7-1.5 km/s and ~+-80-120 ppm; digitisation error ~0.4 km/s, ~30 ppm. -0.2 h to +0.7 h excluded (Doppler shadow).',
                           'depth_definition': 'Gaussian amplitude of CCF residual converted to differential transit depth (ppm)'},
 'stellar_params': ST_W76_E20})
W76_FE.append({
 'id': 'W76_FeCCF_KesseliSnellen2021',
 'citation': {'paper': 'Kesseli & Snellen 2021, ApJL 908, L17', 'arxiv': '2102.01095', 'doi': '10.3847/2041-8213/abe047', 'where': 'Sect. 4 text'},
 'instrument': 'HARPS, 4 transits, Fe I template CCF', 'Kp_kms': 196.52, 'status': 'tabulated',
 'segments': [{'label': 'first half', 'phase_range': [-0.042, -0.026], 'v_kms': V(0.04, 1.47), 'v_kms_alt': V(0.3, 1.07)},
              {'label': 'second half', 'phase_range': [0.002, 0.042], 'v_kms': V(-6.9, 0.75), 'v_kms_alt': V(-7.1, 0.74)}],
 'depth': 'mean residual amplitude ~500 ppm; second half ~100 ppm stronger than first', 'stellar_params': ST_W76_E20})
W76_FE.append({
 'id': 'W76_FeCCF_Kesseli2022',
 'citation': {'paper': 'Kesseli et al. 2022, AJ (atomic spectral survey of WASP-76b)', 'arxiv': '2111.09916', 'doi': '10.3847/1538-3881/ac4336', 'where': 'Table 2'},
 'instrument': 'ESPRESSO (archival E20 transits), per-species template CCF', 'status': 'tabulated',
 'segments': [{'label': 'RV1 (beginning)', 'phase_range': [-0.04, -0.02], 'v_kms': V(-3.89, 0.48), 'amplitude_ppm': V(440, 47)},
              {'label': 'RV2 (end)', 'phase_range': [0.02, 0.04], 'v_kms': V(-10.26, 0.24), 'amplitude_ppm': V(833, 83)}],
 'Kp_best_snr_kms': V(188, ep=10, em=16), 'FWHM_kms': V(9.05, 0.5), 'peak_SNR': 14.24,
 'Kp_for_restframe': 'E20 orbital solution (196.5 km/s) presumably - not explicitly re-stated', 'stellar_params': ST_W76_E20,
 'notes': 'Other species in the same table (V, Cr, Mn, Ni, Co ...) show the same asymmetry.'})
W76_FE.append({
 'id': 'W76_FeRetrieval_Gandhi2022',
 'citation': {'paper': 'Gandhi et al. 2022, MNRAS', 'arxiv': '2206.11268', 'doi': '10.1093/mnras/stac1744', 'where': 'Sect. 5.3 text, Table 2'},
 'instrument': 'ESPRESSO (E20 data), HyDRA-2D retrieval (morning/evening)', 'Kp_kms': 196, 'status': 'tabulated',
 'segments': [{'label': 'first quarter', 'phase_range': [-0.04, -0.02], 'dVsys_kms': V(-5.9, ep=1.5, em=1.1), 'dVwind_FWHM_kms': V(10.2, 1.0),
               'logX_Fe_evening': V(-4.97, ep=0.63, em=0.56), 'logX_Fe_morning': V(-3.83, ep=0.51, em=0.59)},
              {'label': 'last quarter', 'phase_range': [0.02, 0.04], 'dVsys_kms': V(-9.8, ep=1.2, em=1.1),
               'logX_Fe_evening': V(-4.03, ep=0.28, em=0.31), 'logX_Fe_morning': V(-4.59, ep=0.85, em=1.0)}],
 'notes': 'dVsys = retrieved day-night wind shift (one value per phase range, both limbs).'})
W76_FE.append({
 'id': 'W76_FeRetrieval_Gandhi2023',
 'citation': {'paper': 'Gandhi et al. 2023, AJ 165, 242', 'arxiv': '2305.17228', 'doi': '10.3847/1538-3881/accd65', 'where': 'Table 5 (winds), Table 3 (abundances), Table 1 (params)'},
 'instrument': 'ESPRESSO (2 E20 transits), HyDRA-H / HyDRA-2D', 'Kp_kms': V(196.5, 1), 'vsys_kms': -1.16, 'status': 'tabulated',
 'segments': [{'label': 'combined', 'phase_range': [[-0.04, -0.02], [0.02, 0.04]], 'dVsys_kms': V(-6.87, 0.17), 'dVwind_kms': V(6.93, ep=0.55, em=0.5)},
              {'label': '2D first', 'phase_range': [-0.04, -0.02], 'dVsys_kms': V(-4.94, ep=1.62, em=1.55), 'dVwind_kms': V(11.58, ep=0.97, em=0.93)},
              {'label': '2D last', 'phase_range': [0.02, 0.04], 'dVsys_kms': V(-9.03, ep=1.03, em=1.01), 'dVwind_kms': V(6.46, ep=0.55, em=0.56)}],
 'logX_Fe': V(-4.24, ep=0.26, em=0.24), 'stellar_params': {'Rstar_Rsun': 1.756, 'FeH': 0.366}})
W76_FE.append({
 'id': 'W76_Pelletier2023',
 'citation': {'paper': 'Pelletier et al. 2023, Nature 619, 491', 'arxiv': '2306.08739', 'doi': '10.1038/s41586-023-06134-0', 'where': 'Fig. 3, Methods'},
 'instrument': 'MAROON-X (Gemini-N), 3 transits', 'status': 'qualitative only',
 'summary': {'retrieved_Kp_kms': V(180.7, 0.6), 'retrieved_Vsys_kms': V(-7.0, 0.1),
             'trail': 'Fe (and 12 other species) trail starts near 0 km/s and blueshifts over first half of transit, then holds steady (kink); no per-phase table'},
 'notes': 'Kp/Vsys offsets reflect the kinked trail, not true orbital values. Fig. 3 trail not digitised.'})

# ================================================================ Fe I transmission, WASP-121b
s25_even = [(0.0418, -8.77), (0.0386, -6.86), (0.0354, -6.86), (0.0321, -6.53), (0.0289, -6.68), (0.0255, -7.83), (0.0223, -6.86)]
s25_morn = [(-0.0204, -4.55), (-0.0234, -3.86), (-0.0264, -3.29), (-0.0295, -4.33), (-0.0325, -3.72), (-0.0354, -5.88), (-0.0384, -0.72), (-0.0414, -1.62), (-0.0445, 1.77)]
s25 = sorted(s25_morn + s25_even)
W121_FE.append({
 'id': 'W121_FeCCF_Seidel2025',
 'citation': {'paper': 'Seidel et al. 2025, Nature (Vertical structure of an exoplanet atmospheric jet stream)', 'arxiv': '2502.12261',
              'doi': '10.1038/s41586-025-08664-1', 'where': 'Extended Data Table (track_comparison) + Fig. 4a (digitised)'},
 'instrument': 'ESPRESSO 4-UT, two partial transits (morning + evening segments), Fe I template CCF', 'status': 'tabulated segments + digitised per-exposure',
 'segments': [{'label': 'morning segment (pre-mid-transit)', 'phase_range': [-0.0475, -0.017], 'v_kms': V(-4.12, 0.15), 'width_kms': V(15.0, 0.4)},
              {'label': 'evening segment (post-mid-transit)', 'phase_range': [0.017, 0.0475], 'v_kms': V(-6.90, 0.11), 'width_kms': V(13.1, 0.3)}],
 'segments_full_trace_planet_template': [{'label': 'morning', 'v_kms': V(-2.68, 0.14), 'width_kms': V(11.9, 0.3)},
                                         {'label': 'evening', 'v_kms': V(-7.51, 0.16), 'width_kms': V(10.7, 0.3)}],
 'time_series_digitised': {'phase': [a for a, b in s25], 'v_kms': [b for a, b in s25], 'v_err_kms': 0.5,
                           'note': 'dark-purple points in Fig. 4a identified as the Fe I trace by matching segment means (caption gives no colour key); errors ~0.3-0.8 km/s; |phase|<0.017 masked (RM)'},
 'derived_wind_kms': 'Fe I: day-to-night flow ~ -6 to -10 km/s (de-projected); full trace evening: 9.8+-1.6 km/s',
 'depth': 'not tabulated for Fe I',
 'stellar_params': {'Rstar_Rsun': V(1.44, 0.03), 'FeH': V(0.13, 0.04), 'RpRs': V(0.1227, 0.0011), 'vsys_kms': 38.198, 'source': 'Seidel 2025 Table 1 (Borsa 2021 star)'},
 'notes': 'Width = Gaussian width (km/s) as listed (likely sigma or FWHM - table header only says width).'})
W121_FE.append({
 'id': 'W121_CCF_Borsa2021',
 'citation': {'paper': 'Borsa et al. 2021, A&A 645, A24', 'arxiv': '2011.01245', 'doi': '10.1051/0004-6361/202039344', 'where': 'Sect. 5.2 text, Fig. 6'},
 'instrument': 'ESPRESSO 1-UT full transit + 4-UT partial; stellar F6/F9 mask CCF (~76% Fe I weight)', 'status': 'tabulated halves',
 'segments': [{'label': 'first half (excl. ingress)', 'v_kms': V(-2.80, 0.28), 'contrast_ppm': V(1050, 80), 'FWHM_kms': V(7.3, 0.7)},
              {'label': 'second half (excl. egress)', 'v_kms': V(-7.66, 0.16), 'contrast_ppm': V(970, 30), 'FWHM_kms': V(10.2, 0.4)}],
 'average_blueshift_ARM_kms': V(-4, 1), 'R_atmo_over_Rp': V(1.052, 0.015),
 'Kp_kms': 'computed from Table 1 orbit (a/R*=3.8131, P, i=88.49, R*=1.44) ~ 220 km/s (not stated)',
 'stellar_params': {'Teff_K': V(6586, 59), 'logg': V(4.47, 0.08), 'FeH': V(0.13, 0.04), 'Rstar_Rsun': V(1.44, 0.03), 'RpRs': V(0.12534, 0.00005), 'vsys_kms': 38.198}})
W121_FE.append({
 'id': 'W121_FeRetrieval_Gandhi2023',
 'citation': {'paper': 'Gandhi et al. 2023, AJ 165, 242', 'arxiv': '2305.17228', 'doi': '10.3847/1538-3881/accd65', 'where': 'Table 5'},
 'instrument': 'ESPRESSO (1 transit, 2019-01-06), HyDRA-H/2D', 'Kp_kms': V(221, 4.5), 'vsys_kms': 38.198, 'status': 'tabulated',
 'segments': [{'label': 'combined', 'phase_range': [[-0.0475, -0.02], [0.02, 0.0475]], 'dVsys_kms': V(-5.58, ep=0.4, em=0.44), 'dVwind_kms': V(7.92, ep=1.1, em=1.31)},
              {'label': '2D first', 'phase_range': [-0.0475, -0.02], 'dVsys_kms': V(-7.37, ep=1.46, em=1.34), 'dVwind_kms': V(4.02, ep=1.39, em=1.46)},
              {'label': '2D last', 'phase_range': [0.02, 0.0475], 'dVsys_kms': V(-4.33, ep=1.72, em=1.57), 'dVwind_kms': V(4.75, ep=1.52, em=1.64)}],
 'logX_Fe': V(-4.49, ep=0.44, em=0.36), 'stellar_params': {'Rstar_Rsun': 1.44, 'FeH': 0.13}})
W121_FE.append({
 'id': 'W121_Prinoth2025',
 'citation': {'paper': 'Prinoth et al. 2025, A&A (Titanium chemistry of WASP-121b with ESPRESSO 4-UT)', 'arxiv': '2502.12262', 'doi': '10.1051/0004-6361/202452405'},
 'status': 'no Fe numbers', 'notes': 'Same 4-UT data as Seidel 2025; Fe I detected in time-resolved CCF (figure only); time-resolved velocities tabulated/plotted only for Ti I and V I. Stellar Teff 6790+-96 / 6752+-99 (ingress/egress spectra), R* 1.40+-0.05.'})
W121_FE.append({
 'id': 'W121_FeEmission_Hoeijmakers2024',
 'citation': {'paper': 'Hoeijmakers et al. 2024, A&A (Mantis Network IV, Ti cold trap)', 'arxiv': '2210.12847', 'where': 'Table 3'},
 'instrument': 'ESPRESSO, 8 epochs, DAYSIDE EMISSION (not transmission)', 'status': 'tabulated',
 'FeI_emission': {'A_ppm': V(130, 15), 'gamma': V(1.37, 0.45), 'sigma_w_kms': V(4.21, 0.36), 'theta_deg': V(-1.5, 7.6),
                  'theta_deg_east_converted': V(1.5, 7.6), 'v_orb_kms': V(216.24, 0.51), 'v_sys_kms': V(38.72, 0.37)},
 'notes': 'theta = phase-angle offset of peak line emission; model peaks at phase 0.5+theta, so positive = after eclipse (westward) - sign flipped (assumption from eq. 3).'})

out = {
 'meta': {
  'created': '2026-11-02', 'by': 'literature extraction agent (Caltech session)',
  'conventions': {'phase': '0 = mid-transit, 0.5 = mid-eclipse (all converted to [0,1) unless a signed range is given for transit segments)',
                  'flux': 'Fp/Fs in ppm', 'offset': 'hotspot_offset_deg_east: + = peak before mid-eclipse (eastward hotspot); original convention stated per entry',
                  'velocity': 'km/s, planet rest frame, negative = blueshift', 'band_um': 'band edges in micron (Spitzer edges are nominal IRAC ranges)',
                  'values': "{'v': value, 'e': sym 1-sigma} or {'v', 'ep': +err, 'em': -err (magnitude)}", 'halfwidth_um': 'half bin width'},
  'status_flags': 'tabulated = transcribed from table/text or data release; derived = computed here from tabulated data; digitised = read off a figure (approx.)',
  'lit_dir': LIT, 'notes_file': '/resnick/groups/carnegie_poc/jingze/obscheck_1102/obsdata_notes.md'},
 'WASP-121b': {'emission': W121, 'Fe_transmission_hires': W121_FE,
               'system_reference': {'P_d': 1.27492504, 'Kp_kms_Sing2025': V(215.7, 1.1), 'Rp_RJ_Sing2025': V(1.742, 0.006), 'Teq_K_Sing2025': V(2409, 24)}},
 'WASP-76b': {'emission': W76, 'Fe_transmission_hires': W76_FE,
              'system_reference': {'P_d': 1.80988198, 'Kp_kms_E20': V(196.52, 0.94), 'transit_duration_h_approx': 3.75}},
}
json.dump(out, open(OUT, 'w'), indent=1)
print('wrote', OUT, os.path.getsize(OUT))
