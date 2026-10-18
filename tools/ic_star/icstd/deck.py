"""Deck block (athinput fragment) for one standard IC, tools/ic_star.

Emits every star/IC key and the IC-chain/background keys of the AG Car E1 standard
(/raven/ptmp/jinma/overlap_1015/decks/e1_prod_rol2u.athinput md5 552f8769).  Keys NOT emitted
(they are the run's, not the star's): the radial grid (nx1, f_stretch_*), the <rad_m1> solver
set, outputs, tlim (long-tlim rule), vet_gd_rotate_dt (~20 dt).  Paste the blocks into a copy of
the E1 deck (each key replaces the deck's line of the same name in the same block).
"""


def g(x):
    return '%.10e' % x


def block(v):
    """v: dict of resolved values (see make_star_ic.py)"""
    L = []

    def k(blk, key, val, com=''):
        L.append((blk, key, val, com))
    c = v
    k('mesh', 'x1min', g(c['x1min']), '%s: %.3f Rsun' % (c['rin_rule'], c['x1min']/6.957e10))
    k('mesh', 'x1max', g(c['x1max']), '%g R_ph (r_top rule)' % c['rtop_fac'])
    k('mesh', 'ix1_bc', 'user', '')
    k('mesh', 'ox1_bc', 'user', '')
    k('hydro', 'eos', 'general', 'general EOS table (Saha H/He, inert metals)')
    k('hydro', 'general_eos', 'table', '')
    k('hydro', 'eos_h2', 'false', 'H2 OFF (see regime flags if Teff < 5000 K)')
    k('hydro', 'eos_ionization', 'true', '')
    k('hydro', 'eos_metal_ionization', 'false', '')
    k('hydro', 'eos_radiation', 'false', 'MANDATORY with <rad_m1>')
    k('hydro', 'eos_xh', repr(c['X']), '')
    k('hydro', 'eos_yhe', repr(c['Y']), 'Z = %g' % c['Z'])
    k('hydro', 'eos_a_metal', repr(c['a_metal']), 'inert-metal mean weight (E1 default 16)')
    k('hydro', 'eos_logd_min', '-21.0', '= he_opac_logd_min, covers dfloor')
    k('hydro', 'eos_logd_max', repr(c['eos_logd_max']), '')
    k('hydro', 'eos_logt_min', '%.1f' % c['eos_logt_min'], 'floor(10 log10 tfloor)/10')
    k('hydro', 'eos_logt_max', repr(c['eos_logt_max']), 'column max T %.3e K' % c['Tmax'])
    k('hydro', 'eos_dlog', repr(c['eos_dlog']), '')
    k('hydro', 'wb_rmax', g(c['wb_rmax']), '%g R_ph' % c['wb_rmax_rel'])
    k('hydro', 'dfloor', '%.1e' % c['dfloor'], '%g x column floor' % c['dfloor_rel'])
    k('hydro', 'tfloor_kelvin', '%.1f' % c['tfloor_kelvin'], c['why_tfloor'])
    k('hydro', 'pfloor', '%.1e' % c['pfloor'], c['why_pfloor'])
    k('hydro', 'dfloor_keep_velocity', 'true', '')
    k('hydro', 'vceil', '%.1e' % c['vceil'], c['why_vceil'])
    k('rad_m1', 'implicit_bc_x1min', 'flux', '')
    k('rad_m1', 'implicit_flux_x1min', '%.12e' % c['F_rin'],
      'column F_r(x1min), log-linear (he_star_m1 check 1e-3); L/(4 pi x1min^2) = %.12e'
      % c['F_L_rin'])
    k('rad_m1', 'implicit_bc_x1max', 'marshak', '')
    k('rad_m1', 'implicit_flux_x1max', '0.0', '')
    p = 'problem'
    k(p, 'he_gm', g(c['GM']), c['gm_src'])
    if c.get('gm_column'):
        k(p, 'he_gm_column', 'true', "mode mesa: G m(r) = he_gm + G int 4 pi r^2 rho_IC dr")
    else:
        k(p, 'he_gm_column', 'false', 'mode global: point mass, as the column')
    k(p, 'he_ic_file', c['ic_file'], 'md5 %s' % c['ic_md5'])
    k(p, 'he_opac_table', c['ross'], 'md5 %s' % c['ross_md5'])
    k(p, 'he_planck_table', c['planck'], 'md5 %s' % c['planck_md5'])
    k(p, 'he_opac_logd_min', '-21.0', '')
    k(p, 'he_nfine', '262144', '')
    k(p, 'he_tau_int', '1.0', '')
    k(p, 'he_ic_cols', '5', '7-column file')
    k(p, 'he_ic_rad', 'sph1d', 'E1 IC chain')
    k(p, 'he_ic_rad_gas', 're_atm', '')
    k(p, 'he_ic_rad_hse', 'true', '')
    k(p, 'he_ic_rad_atm_force', 'module', '')
    k(p, 'he_ic_rad_atm_post', 'false', 'exclusive with he_ic_equil')
    k(p, 'he_ic_eint_from_t', 'true', 'eint = e_code(rho, T col 6)')
    k(p, 'he_ic_balance', 'true', '')
    k(p, 'he_ic_balance_rmax', g(c['he_ic_balance_rmax']),
      '%g R_ph (inside the MLT++ adiabat, below the tau-3 join %.4f R_ph)'
      % (c['balance_rmax_rel'], c['rj']/c['R']))
    k(p, 'he_ic_equil', 'true', '')
    k(p, 'he_ic_equil_dtfac', '1.0', '')
    k(p, 'he_ic_equil_relax_n', '40', '')
    k(p, 'he_ic_equil_relax_tol', '1.0e-4', '')
    k(p, 'he_ic_equil_maxit', '1', '')
    k(p, 'he_ic_equil_tol', '1.0e-9', '')
    k(p, 'he_ic_equil_last_march', 'true', '')
    k(p, 'he_wall_zero_flux', 'true', '')
    k(p, 'he_wall_closed', 'true', '')
    k(p, 'he_bc_inner', 'wall', '')
    k(p, 'mlt_flux_frozen', 'true', 'frozen MLT++ scaffold (col 7)')
    k(p, 'mlt_closure', 'frozen', '')
    if c.get('t_turn') and c['t_turn'] == c['t_turn']:
        k(p, 'mlt_ramp_start', '%.4e' % c['ramp_start'],
          '%g x t_Fe (H_p/v_MLT at the Fe peak = %.4e s)' % (c['start_turnovers'], c['t_turn']))
        k(p, 'mlt_ramp_time', '%.4e' % c['ramp_time'], '%g x t_Fe' % c['time_turnovers'])
    else:
        k(p, 'mlt_ramp_start', '-1.0', 'no convective Fe zone in the column: no ramp')
        k(p, 'mlt_ramp_time', '0.0', '')
    k(p, 'he_bc_outer', 'outflow', 'outflow top, no inflow')
    k(p, 'he_sponge_rate', '%.1e' % c['he_sponge_rate'], c['why_sponge'])
    k(p, 'he_sponge_r0', g(c['he_sponge_r0']), '%g x1max' % c['sponge_r0_rel'])
    k(p, 'he_sponge_mode', 'radial', '')
    k(p, 'he_sponge_dmax', '%.1e' % c['he_sponge_dmax'],
      '%g x column floor %.1e (sponge ACTIVE on floor gas)' % (c['dmax_rel'], c['rho_floor']))
    k(p, 'he_seed_rlo', g(c['seed_rlo']), 'CZ base x 1.007 (E1 146.0/145.0 Rsun)')
    k(p, 'he_seed_rhi', g(c['seed_rhi']), '0.9786 R_ph (E1 380/388.3)')
    k(p, 'he_grid_dump', 'true', '')
    for kv in c.get('deck_extra', []):      # spec deck_extra: 'problem/key=value' (opt-in keys)
        b_, kv_ = kv.split('/', 1)
        k(b_, kv_.split('=')[0], kv_.split('=', 1)[1], 'spec deck_extra')
    out = ['# ---- tools/ic_star deck block for %s (spec %s)' % (c['name'], c['spec']),
           '# ---- star/IC keys of the standard; grid (nx1, f_stretch_*), <rad_m1> solver set,'
           ' outputs, tlim, he_seed amplitude/k band from the E1 deck']
    cur = None
    for blk, key, val, com in L:
        if blk != cur:
            out.append('<%s>' % blk)
            cur = blk
        out.append('%-28s = %-24s%s' % (key, val, ('   # ' + com) if com else ''))
    return '\n'.join(out) + '\n'
