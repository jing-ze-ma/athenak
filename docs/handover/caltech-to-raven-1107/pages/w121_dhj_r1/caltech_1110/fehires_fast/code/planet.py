"""planet.py: per-planet system / instrument / composition parameters for fe_hires.

  cfg = planet.load('planets/wasp121b.json')
Every key has the WASP-121b value (the 10-07 viper result, FE_HIRES_w121_results.md) as its
default, so a JSON only needs the keys that differ.  Derived quantities are added:
  RS, RP [cm], PROT [s], OMEGA [1/s], KP0 [km/s], PH14 (T1/T4), PH23 (T2/T3), LSF [km/s].
"""
import json
import os
import numpy as np

RJ, RSUN, MJ = 7.1492e9, 6.957e10, 1.898e30
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULTS = {
    'name': 'WASP-121b',
    # system
    'R_star_Rsun': 1.461,          # R*
    'R_p_RJ': 1.742,               # Rp (sets the transit contacts and the static test radius)
    'M_p_MJ': 1.17,                # used only by the static=1 validation atmosphere
    'a_over_Rs': 3.7844,           # a/R*
    'b': 0.160,                    # impact parameter
    'P_orb_d': 1.27492504,         # orbital period (= rotation period, tidally locked)
    'omega': 5.704026e-05,         # planet rotation rate [1/s]; null -> 2 pi / P_orb
    'Kp_kms': None,                # null -> 2 pi a/P sqrt(1 - (b/a)^2) (= 219.2 for WASP-121b)
    'limb_darkening': [0.0, 0.0],  # quadratic u1, u2 of the stellar disc; [0, 0] = uniform
    # instrument
    'lsf_fwhm_kms': 4.28,          # Gaussian LSF FWHM (ESPRESSO 4-UT MR); or set 'resolution_R'
    'resolution_R': None,          # if set (and lsf_fwhm_kms null): FWHM = c/R
    # exposures
    'nph': 31,                     # equally spaced exposures
    'phase_range': None,           # [phi0, phi1]; null -> T1..T4 (+-PH14)
    'rm_mask_phase': 0.017,        # |phi| < this removed in the 'rmcut' Kp-Vsys set and segments
    # composition (must match the run's EOS keys; comptab built by comptab.cpp)
    'eos_xh': 0.7188,
    'eos_yhe': 0.2420,
    'comptab': 'data/comptab_3x.bin.gz',
    'a_fe': 3.2e-5,                # Fe/H by number (EOS donor table value; not scaled by [M/H])
    'pcut_bar': 1e-6,              # Fe removed at p < pcut (model top); 0 = no cut
    # ana.py Kp-marg (Seidel+2026 marginalise over the orbital velocity; no range given)
    'kp_marg_range': None,         # [Kp_lo, Kp_hi] km/s, flat prior; null -> Kp0 +- 50
    'kp_marg_snr': None,           # null: noise = model map std at |v| > 30; else map max = this S/N
}


def load(fn=None):
    cfg = dict(DEFAULTS)
    if fn:
        with open(fn) as f:
            user = json.load(f)
        unknown = set(user) - set(DEFAULTS) - {'comment', 'refs'}
        if unknown:
            raise SystemExit('unknown planet keys: %s' % sorted(unknown))
        cfg.update(user)
    cfg['RS'] = cfg['R_star_Rsun']*RSUN
    cfg['RP'] = cfg['R_p_RJ']*RJ
    cfg['MP'] = cfg['M_p_MJ']*MJ
    cfg['PROT'] = cfg['P_orb_d']*86400.0
    cfg['OMEGA'] = cfg['omega'] if cfg['omega'] else 2*np.pi/cfg['PROT']
    aor, b = cfg['a_over_Rs'], cfg['b']
    kp = 2*np.pi*aor*cfg['RS']/cfg['PROT']*np.sqrt(1 - (b/aor)**2)/1e5
    cfg['KP0'] = cfg['Kp_kms'] if cfg['Kp_kms'] else kp
    krp = cfg['RP']/cfg['RS']
    cfg['KRP'] = krp
    cfg['PH14'] = np.arcsin(np.sqrt((1 + krp)**2 - b**2)/aor)/(2*np.pi)
    cfg['PH23'] = np.arcsin(np.sqrt((1 - krp)**2 - b**2)/aor)/(2*np.pi)
    if cfg['lsf_fwhm_kms']:
        cfg['LSF'] = float(cfg['lsf_fwhm_kms'])
    else:
        cfg['LSF'] = 2.99792458e5/float(cfg['resolution_R'])
    if not os.path.isabs(cfg['comptab']):
        cfg['comptab'] = os.path.join(HERE, cfg['comptab'])
    return cfg
