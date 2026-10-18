"""Star spec (YAML or JSON) -> resolved inputs, with the standard defaults (tools/ic_star).

Every default below is the AG Car E1 production value (ic_agcar_A_ge_deep70.txt +
e1_prod_rol2u.athinput) or a rule in units of the star that reproduces it.
"""
import copy
import json
import os

import numpy as np

from .micro import GMSUN, G_CODE, LSUN, RSUN, SIG

DEFAULTS = {
    'composition': {'mixture': 'GS98', 'a_metal': 16.0},
    # MLT: MESA r15140 defaults (star/defaults/controls.defaults): mixing_length_alpha 2
    # (:2018), MLT_option 'Cox' (:2191), alt_scale_height_flag .true. (:2148),
    # okay_to_reduce_gradT_excess .false. (:2292, i.e. MLT++ OFF), gradT_excess_f1 1e-4,
    # f2 1e-3 (:2304-2305), lambda1 1.0 beta1 0.35 (:2342-2343), lambda2 0.5 beta2 0.25
    # (:2355-2356), dlambda 0.1 dbeta 0.1 (:2368-2369), max_logT 8 (:2399).
    # mltpp 'mesa_gated': MESA's (beta, lambda) gate, alpha* = rule(column(alpha*));
    # 'agcar_legacy': a FIXED pp_alpha applied wherever nabla - nabla_ad > f1 (pp_form 'mesa'
    # = MESA's af blend, the AG Car A builder; 'f2' = make_ic_agcar_mltpp.py).
    'mlt': {'alpha': 2.0, 'mlt_option': 'Cox', 'alt_scale_height': True,
            'mltpp': 'mesa_gated', 'use_mltpp': False, 'f1': 1.0e-4, 'f2': 1.0e-3,
            'lambda1': 1.0, 'beta1': 0.35, 'lambda2': 0.5, 'beta2': 0.25, 'dlambda': 0.1,
            'dbeta': 0.1, 'xa1_center': 0.0, 'max_logT': 8.0, 'pp_alpha': 1.0,
            'pp_form': 'mesa', 'pp_tol': 2.0e-3},
    'grid': {'dlnr': 2.0e-4, 'drn_atm_rph': 2.0e-4, 'fl_step': 50.0, 'tauj': 3.0,
             'rtol': 1.0e-12, 'rho_lo': 1.0e-14, 'rho_hi': 1.0e-7, 'atm_rad': 0,
             'scheme': 'agcar_ge'},
    # r_in: the WALL radius x1min.  rule 'hp_below_cz' (default, n_hp 4.5) or 'explicit'
    # (r_rsun or r_cm).  Column bottom = rlo_fac x1min unless bottom_rsun is given.
    'r_in': {'rule': 'hp_below_cz', 'n_hp': 4.5, 'rlo_fac': 0.95, 'cz_fmlt': 0.01,
             'explore_tmax': 3.0e6, 'l_flat_tol': 0.01},
    # r_top: x1max = rtop_fac R_ph; column top = col_top_fac x1max
    'r_top': {'rtop_fac': 3.0, 'col_top_fac': 1.15},
    # background / floor layer rules (hot-star calibrated: they reproduce the AG Car E1 values)
    'background': {
        'rho_floor_rule': 'rel_rhoph', 'rho_floor_rel': 4.0e-6,     # x rho_ph, -> decade
        'dmax_rel': 2.0,                # he_sponge_dmax = dmax_rel x column floor
        'dfloor_rel': 1.0e-4,           # dfloor = dfloor_rel x column floor
        'sponge_rate_rel': 80.0,        # he_sponge_rate = rel x (g/c_s,iso)_ph -> decade
        'sponge_r0_rel': 0.9,           # x x1max
        'tfloor_rel': 0.76,             # tfloor = rel x T_rad(x1max), floored to 100 K
        'wb_rmax_rel': 1.1,             # x R_ph
        'pfloor_rel': 6.0e-8,           # x P_gas(column floor, tfloor) -> decade
        'vceil_rel': 4.3,               # x v_esc(R_ph) -> decade
        'balance_rmax_rel': 0.8808,     # he_ic_balance_rmax x R_ph (E1 342/388.3 Rsun)
    },
    'ramp': {'start_turnovers': 2.0, 'time_turnovers': 10.0},
    # 10-16 seam removal (icstd/seam.py): F_MLT taper tau 30 -> 3 (on for new stars) and the
    # code-consistent atmosphere through the code's own sph1d (needs seam.deck = the run deck)
    'taper': {'enable': True, 'tau_lo': 3.0, 'tau_hi': 30.0},
    'seam': {'enable': True, 'deck': None,
             'binary': '/raven/ptmp/jinma/icstd_1016/bin/athena_he_cpu_b16b66db',
             'tol_cells': 0.1, 'maxit': 8, 'omp': 36, 'timeout': 7200,
             'launcher': ['bash', '-lc', 'module purge >/dev/null 2>&1; '
                          'module load gcc/13 openmpi/5.0; exec "$0" "$@"']},
    'eos': {'logd_min': -21.0, 'logd_max': -1.0, 'logt_max': 6.8, 'dlog': 0.02},
    'opacity': {'ld_ext': -21.0},
}

CFG = {
    'cache': '/raven/ptmp/jinma/ic_star_cache',
    'eos_binary': '/raven/ptmp/jinma/bubbles_1012/bin/athena_cpu_m4',
}


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load(fn):
    txt = open(fn).read()
    if fn.endswith('.json'):
        s = json.loads(txt)
    else:
        try:
            import yaml
        except ImportError as e:
            raise SystemExit('PyYAML missing: use a .json spec') from e
        s = yaml.safe_load(txt)
    s = _merge(DEFAULTS, s)
    s['_cfg'] = _merge(CFG, s.get('config', {}))
    s['_file'] = os.path.abspath(fn)
    return s


def rsun_or_cm(d, key):
    """d[key + '_rsun'] or d[key + '_cm'] -> cm (None if neither)"""
    if d.get(key + '_cm') is not None:
        return float(d[key + '_cm'])
    if d.get(key + '_rsun') is not None:
        return float(d[key + '_rsun'])*RSUN
    return None


def resolve_star(s, mesa=None):
    """L, GM, R_ph, Teff from the spec (and the MESA profile where the spec says 'mesa')."""
    st = s['star']
    out, src = {}, {}
    # L
    if st.get('L') == 'mesa':
        out['L'], src['L'] = mesa['L_ph'], 'MESA photosphere_L'
    elif st.get('L') is not None:
        out['L'], src['L'] = float(st['L']), 'spec L'
    else:
        out['L'], src['L'] = LSUN*10.0**float(st['logL']), 'spec logL (Lsun 3.828e33)'
    # GM
    if st.get('GM') is not None:
        out['GM'], src['GM'] = float(st['GM']), 'spec GM'
    elif st.get('M') == 'mesa':
        out['GM'], src['GM'] = G_CODE*mesa['M_tot'], 'G_code x MESA total mass'
    else:
        out['GM'], src['GM'] = float(st['M'])*GMSUN, 'spec M x GMsun 1.3271244e26'
    # R_ph / Teff
    if st.get('R_ph') == 'mesa':
        out['R'], src['R'] = mesa['R_ph'], 'MESA photosphere_r'
    elif st.get('R_ph_cm') is not None or st.get('R_ph_rsun') is not None:
        out['R'], src['R'] = rsun_or_cm(st, 'R_ph'), 'spec R_ph'
    else:
        teff = mesa['Teff'] if st.get('Teff') == 'mesa' else float(st['Teff'])
        out['R'] = np.sqrt(out['L']/(4*np.pi*SIG*teff**4))
        src['R'] = 'sqrt(L/(4 pi sigma Teff^4)), Teff %s' % teff
    out['Teff'] = (out['L']/(4*np.pi*out['R']**2*SIG))**0.25
    return out, src
