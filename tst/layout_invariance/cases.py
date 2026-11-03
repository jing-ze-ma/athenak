"""cases.py: the layout-invariance cases run by li.py.

Each case: name, bin (family in bins.cfg), deck (in decks/), set (deck overrides, block/key),
vars (output variables, one bin file each), ncyc, layouts (the FIRST is the reference; it is
also run with the 1-ulp perturbation), expect (PASS | FAIL: a known, documented layout
dependence that the case keeps visible), what (one line).  Layout: name, mb=(nx1,nx2,nx3)
of the MeshBlock, np (MPI ranks), set (extra overrides for that layout only).
Keep the grids tiny: the whole suite is meant to run in minutes on <= 16 CPU cores."""
import os

PKG = os.environ.get('LI_DHJ_PKG', '/resnick/groups/carnegie_poc/jingze/w121r1_1101/pkg')


def L(name, mb, np=1, **kw):
    d = dict(name=name, mb=mb, np=np)
    d.update(kw)
    return d


def cs_layouts(n, nr, rad_split=True):
    """n = cells per panel side, nr = radial cells: 1 block/panel (ref), 2x2, 4x4 per panel,
    radial split, and rank counts that cut panels between ranks."""
    h, q = n // 2, n // 4
    lays = [L('b1x1_np1', (nr, n, n)),
            L('b2x2_np1', (nr, h, h)),
            L('b2x2_np4', (nr, h, h), 4),
            L('b4x4_np6', (nr, q, q), 6),
            L('b1x1_np3', (nr, n, n), 3)]
    if rad_split:
        lays.append(L('b2x2r2_np2', (nr // 2, h, h), 2))
    return lays


def sp_layouts(nr, nth, nph, rsplit=True):
    """spherical polar: the reference has 2 blocks in phi (the polar boundary needs an
    even number); then 4 and 8 in phi (the production accretor has 8), theta and radial
    splits, several ranks"""
    return [L('ph2_np1', (nr, nth, nph // 2)),
            L('ph4_np1', (nr, nth, nph // 4)),
            L('ph8_np4', (nr, nth, nph // 8), 4),
            L('th2ph2_np2', (nr, nth // 2, nph // 2), 2),
            L('r2ph4_np8', (nr // 2, nth, nph // 4), 8)][:5 if rsplit else 4]


def cart_layouts(n1, n2, n3):
    lays = [L('b1_np1', (n1, n2, n3)),
            L('b2_np1', (n1 // 2, n2, n3)),
            L('b8_np1', (n1 // 2, max(n2 // 2, 1), max(n3 // 2, 1)))]
    lays.append(L('b8_np4', lays[-1]['mb'], 4))
    if n3 > 1:
        lays.append(L('b16_np8', (n1 // 4, n2 // 2, n3 // 2), 8))
    return lays


CS_HYD = dict(bin='cs', deck='cubed_sphere_blast.athinput', ncyc=10,
              vars=['hydro_u'],
              set={'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16, 'hydro/nscalars': 1})
CS_MHD = dict(bin='cs', deck='cubed_sphere_mhd_blast.athinput', ncyc=10,
              vars=['mhd_u', 'mhd_bcc'],
              set={'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16, 'mhd/nscalars': 1})


def var(base, name, what, expect='PASS', extra=None, **kw):
    c = dict(base)
    c['set'] = dict(base['set'])
    c['set'].update(extra or {})
    c.update(name=name, what=what, expect=expect)
    c.update(kw)
    if 'layouts' not in c:
        n = int(c['set'].get('mesh/nx2', 16))
        c['layouts'] = cs_layouts(n, int(c['set'].get('mesh/nx1', 8)))
    return c


CASES = [
    # ------------------------------------------------------------ Cartesian (built-in pgens)
    dict(name='cart_lw_hydro', bin='bip', deck='cart_lw_hydro.athinput', ncyc=10,
         vars=['hydro_u'], what='3-D hydro linear wave, periodic, plm/llf',
         set={'mesh/nx1': 32, 'mesh/nx2': 16, 'mesh/nx3': 16, 'hydro/rsolver': 'hllc'},
         layouts=cart_layouts(32, 16, 16)),
    dict(name='cart_implode', bin='bip', deck='cart_implode.athinput', ncyc=20,
         vars=['hydro_u'], what='2-D implode, REFLECTING walls (physical-BC corners)',
         set={'mesh/nx1': 32, 'mesh/nx2': 32}, layouts=cart_layouts(32, 32, 1)),
    dict(name='cart_lw_mhd', bin='bip', deck='cart_lw_mhd.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'], what='3-D MHD linear wave, CT edges shared by 4 blocks',
         set={'mesh/nx1': 32, 'mesh/nx2': 16, 'mesh/nx3': 16, 'mhd/rsolver': 'hlld'},
         layouts=cart_layouts(32, 16, 16)),
    dict(name='cart_ot', bin='bip', deck='cart_ot.athinput', ncyc=20,
         vars=['mhd_u', 'mhd_bcc'], what='2-D Orszag-Tang MHD',
         set={'mesh/nx1': 32, 'mesh/nx2': 32}, layouts=cart_layouts(32, 32, 1)),
    dict(name='cart_m1_beam', bin='bip', deck='cart_m1_beam.athinput', ncyc=10,
         vars=['m1'], what='2-D M1 beam, user/vacuum boundaries',
         set={'mesh/nx1': 32, 'mesh/nx2': 32}, layouts=cart_layouts(32, 32, 1),
         perturb='rad_m1/cfl_rad'),
    # ------------------------------------------------------------ cubed sphere (cs_test)
    var(CS_HYD, 'cs_hyd_clamp', 'cs hydro blast through seam + cube vertex, '
        'cs_seam_blockend = clamp (the 11-03 bug, kept visible)', expect='FAIL'),
    var(CS_HYD, 'cs_hyd_extend', 'cs hydro blast, cs_seam_blockend = extend',
        extra={'mesh/cs_seam_blockend': 'extend'}),
    var(CS_HYD, 'cs_hyd_ng3', 'cs hydro blast, nghost 3 ppm (HG/RSG halo depth), extend',
        extra={'mesh/cs_seam_blockend': 'extend', 'mesh/nghost': 3,
               'hydro/reconstruct': 'ppm4'}),
    var(CS_MHD, 'cs_mhd_clamp', 'cs MHD blast (hlld, CT), clamp', expect='FAIL'),
    var(CS_MHD, 'cs_mhd_extend', 'cs MHD blast (hlld, CT), extend',
        extra={'mesh/cs_seam_blockend': 'extend'}),
    var(CS_MHD, 'cs_mhd_shared_o', 'cs MHD blast, r1m seam keys (cs_seam_flux shared, '
        'cs_seam_mhd shared_o, uct_hlld), extend',
        extra={'mesh/cs_seam_blockend': 'extend', 'mesh/cs_seam_flux': 'shared',
               'mhd/cs_seam_mhd': 'shared_o', 'mhd/emf_scheme': 'uct_hlld'}),
    dict(name='cs_resist_rkg', bin='cs', deck='cubed_sphere_resist.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'], what='cs Ohmic + ambipolar (constant), RKG STS, extend',
         set={'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16,
              'mesh/cs_seam_blockend': 'extend', 'mhd/ambipolar': 'constant',
              'mhd/eta_ad': 0.05, 'mhd/use_rkg_sts': 'true'},
         layouts=cs_layouts(16, 8)),
    dict(name='cs_resist_rkg_clamp', bin='cs', deck='cubed_sphere_resist.athinput',
         ncyc=10, vars=['mhd_u', 'mhd_bcc'], what='as cs_resist_rkg, clamp', expect='FAIL',
         set={'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16, 'mhd/ambipolar': 'constant',
              'mhd/eta_ad': 0.05, 'mhd/use_rkg_sts': 'true'},
         layouts=cs_layouts(16, 8)),
    dict(name='cs_raddiff', bin='cs', deck='cubed_sphere_raddiff.athinput', ncyc=10,
         vars=['hydro_u'], what='cs scalar diffusion (iprob 15), extend',
         set={'mesh/cs_seam_blockend': 'extend'}, layouts=cs_layouts(16, 8)),
    # ------------------------------------------------------------ spherical polar (sp_test)
    dict(name='sp_hyd_blast', bin='sp', deck='spherical_polar_blast.athinput', ncyc=10,
         vars=['hydro_u'], what='sp hydro blast, polar boundary, phi split 2/4/8',
         set={}, layouts=sp_layouts(8, 16, 32)),
    dict(name='sp_hyd_ringavg', bin='sp', deck='spherical_polar_blast.athinput', ncyc=10,
         vars=['hydro_u'], what='sp hydro blast + polar_ring_avg (auto group sizes)',
         set={'mesh/polar_ring_avg': 'true', 'mesh/nx3': 64}, layouts=sp_layouts(8, 16, 64)),
    dict(name='sp_hyd_ringfourier', bin='sp', deck='spherical_polar_blast.athinput',
         ncyc=10, vars=['hydro_u'], what='sp hydro blast + ring Fourier filter',
         set={'mesh/polar_ring_avg': 'true', 'mesh/polar_ring_avg_filter': 'fourier', 'mesh/nx3': 64},
         layouts=sp_layouts(8, 16, 64)),
    dict(name='sp_mhd_blast', bin='sp', deck='spherical_polar_blast_mhd.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'], what='sp MHD blast, polar boundary, phi split',
         set={}, layouts=sp_layouts(8, 16, 32, False)),
    dict(name='sp_mhd_polar_keys', bin='sp', deck='spherical_polar_blast_mhd.athinput',
         ncyc=10, vars=['mhd_u', 'mhd_bcc'],
         what='sp MHD blast + polar_emf_diss + polar_hlle_rows + polar_ring_avg',
         set={'mhd/polar_emf_diss': 'true', 'mhd/polar_hlle_rows': '2',
              'mesh/polar_ring_avg': 'true', 'mesh/nx3': 64},
         layouts=sp_layouts(8, 16, 64, False)),
    dict(name='sp_resist', bin='sp', deck='spherical_polar_resist.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'], what='sp Ohmic + ambipolar decay, RKG, polar',
         set={'mhd/ambipolar': 'constant', 'mhd/eta_ad': 0.05, 'mhd/use_rkg_sts': 'true'},
         layouts=sp_layouts(8, 16, 32, False)),
]

# ---------------------------------------------------------------- deep hot Jupiter (dhj)
# The production decks (r1m v6 MHD, hydro planet) on a C8 x 32 grid.  Needs the W121
# packages (ck tables, IC profile) linked into the run dir: LI_DHJ_PKG.
DHJ_LINKS = {os.path.join(PKG, 'ck'): 'ck', os.path.join(PKG, 'ic'): 'ic'}


def dhj_layouts():
    return [L('b1x1_np1', (32, 8, 8)), L('b2x2_np1', (32, 4, 4)), L('b2x2_np4', (32, 4, 4), 4),
            L('b1x1_np3', (32, 8, 8), 3), L('b1x1_np6', (32, 8, 8), 6)]
    # (no radial split: deep_hot_jupiter_rt refuses it)


DHJ_SET = {'mesh/nx1': 32, 'mesh/nx2': 8, 'mesh/nx3': 8, 'time/ndiag': 1}
DHJ_M = dict(bin='dhj', deck='dhj_r1m.athinput', ncyc=12, vars=['mhd_u', 'mhd_bcc'],
             links=DHJ_LINKS, timeout=7200, set=dict(DHJ_SET), layouts=dhj_layouts())
DHJ_H = dict(bin='dhj', deck='dhj_hyd.athinput', ncyc=12, vars=['hydro_u'],
             links=DHJ_LINKS, timeout=7200, set=dict(DHJ_SET), layouts=dhj_layouts())


def dvar(base, name, what, extra=None, expect='PASS'):
    c = dict(base)
    c['set'] = dict(base['set'])
    c['set'].update(extra or {})
    c.update(name=name, what=what, expect=expect)
    return c


EXT = {'mesh/cs_seam_blockend': 'extend'}
CASES += [
    dvar(DHJ_M, 'dhj_r1m', 'r1m v6 MHD deck (wb, sponge, Roche + wb_tan, outflow1, ck every 4, '
         'RKG Ohmic+ambipolar, tracers), extend', EXT),
    dvar(DHJ_M, 'dhj_r1m_clamp', 'r1m v6 deck, clamp (the prod default)', expect='FAIL'),
    dvar(DHJ_H, 'dhj_hyd', 'hydro planet deck (lhllc), extend', EXT),
    dvar(DHJ_H, 'dhj_hyd_hse_spb', 'hydro deck + top_bc hse + sponge_bottom, extend',
         dict(EXT, **{'problem/top_bc': 'hse', 'problem/sponge_bottom': 'true'})),
    dvar(DHJ_H, 'dhj_hyd_nojreuse', 'hydro deck, ck_impl_jreuse = 0 (control for the per-rank '
         'chord-Jacobian rebuild decision)', dict(EXT, **{'problem/ck_impl_jreuse': 0.0})),
    dvar(DHJ_H, 'dhj_hyd_ck1', 'hydro deck, ck every cycle, xstep 1, extend',
         dict(EXT, **{'problem/ck_impl_every': 1, 'problem/ck_impl_xstep': 1})),
]

# ---------------------------------------------------------------- probes of the extend pass
# (cs_mhd_extend fails for 4x4 blocks per panel at 4-cell blocks): which matters, the number
# of blocks per panel side or the block size (= 2 nghost)?
CASES += [
    var(CS_MHD, 'cs_mhd_ext_n32', 'cs MHD blast C32: 1 block vs 4x4 blocks of 8 cells, extend',
        extra={'mesh/cs_seam_blockend': 'extend', 'mesh/nx2': 32, 'mesh/nx3': 32},
        layouts=[L('b1x1_np1', (8, 32, 32)), L('b4x4_np1', (8, 8, 8)),
                 L('b2x2_np1', (8, 16, 16)), L('b8x8_np6', (8, 4, 4), 6)]),
    var(CS_MHD, 'cs_mhd_ext_n8', 'cs MHD blast C8: 1 block vs 2x2 blocks of 4 cells, extend',
        extra={'mesh/cs_seam_blockend': 'extend', 'mesh/nx2': 8, 'mesh/nx3': 8},
        layouts=[L('b1x1_np1', (8, 8, 8)), L('b2x2_np1', (8, 4, 4))]),
    var(CS_HYD, 'cs_hyd_ext_n32', 'cs hydro blast C32: 1 vs 4x4 / 8x8 blocks, extend',
        extra={'mesh/cs_seam_blockend': 'extend', 'mesh/nx2': 32, 'mesh/nx3': 32},
        layouts=[L('b1x1_np1', (8, 32, 32)), L('b4x4_np1', (8, 8, 8)),
                 L('b8x8_np6', (8, 4, 4), 6)]),
]

# ---------------------------------------------------------------- HG/RSG and the sp accretor
# Production decks (absolute table paths) on tiny grids.  HG (he_star_m1, cubed sphere,
# nghost 3, M1 vet_col + vet_gd twin + vet_sc): C16 x 64, 1 vs 2x2 blocks per panel.
# Accretor (ry_per_accretor, spherical polar, polar ring average _ke/_gz, vet_gd, M1
# hesdirk2): 50 x 29 x 64, phi split 2 / 4 / 8 (production: 8 blocks of 64 in phi).
CASES += [
    dict(name='he_hg', bin='he', deck='he_hg.athinput', ncyc=3, timeout=7200,
         perturb='problem/he_gm',   # dt is set by vet_gd_rotate_dt alignment, not the cfl
         vars=['hydro_u', 'm1'], what='HG production deck (he_star_m1, M1 vet), C16 x 704 (production radial grid)',
         set={'mesh/nx2': 16, 'mesh/nx3': 16, 'mesh/cs_seam_blockend': 'extend'},
         layouts=[L('b1x1_np6', (704, 16, 16), 6), L('b2x2_np6', (704, 8, 8), 6),
                  L('b1x1_np3', (704, 16, 16), 3)]),
    dict(name='acc_sp', bin='acc', deck='acc_sp.athinput', ncyc=3, timeout=7200,
         vars=['hydro_u', 'm1'], what='sp accretor production deck (ring avg ke/gz, vet_gd), '
         '500 x 48 x 64 (plateau core 16), phi split',
         set={'mesh/nx2': 48, 'mesh/nx3': 64, 'mesh/theta_plateau_ncore': 16,
              'mesh/theta_plateau_dth': 0.06},
         layouts=[L('ph2_np2', (500, 48, 32), 2), L('ph4_np4', (500, 48, 16), 4)]),
]
# (8 phi blocks of 8 cells: refused since 8995b8d8e, the pole group is 16)
_acc = [c for c in CASES if c['name'] == 'acc_sp'][0]
CASES.append(dict(_acc, name='acc_sp_novetgd', set=dict(_acc['set'], **{'rad_m1/vet_gd': 'false'}),
                  what='accretor deck with rad_m1/vet_gd = false (control: is the M1 layout '
                  'dependence vet_gd?)'))

# ---------------------------------------------------------------- resistive probes
# (sp_resist and cs_resist_rkg fail with any block split): which term?
_SPR = dict(bin='sp', deck='spherical_polar_resist.athinput', ncyc=10, vars=['mhd_u', 'mhd_bcc'],
            layouts=sp_layouts(8, 16, 32, False)[:3])
_CSR = dict(bin='cs', deck='cubed_sphere_resist.athinput', ncyc=10, vars=['mhd_u', 'mhd_bcc'],
            layouts=cs_layouts(16, 8)[:2])
_CSRS = {'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16, 'mesh/cs_seam_blockend': 'extend'}
for _nm, _s in (('ohm', {}), ('ohm_rkg', {'mhd/use_rkg_sts': 'true'}),
                ('amb', {'mhd/ambipolar': 'constant', 'mhd/eta_ad': 0.05,
                         'mhd/ohmic_resistivity': 'none'}),
                ('amb_rkg', {'mhd/ambipolar': 'constant', 'mhd/eta_ad': 0.05,
                             'mhd/ohmic_resistivity': 'none', 'mhd/use_rkg_sts': 'true'})):
    CASES.append(dict(_SPR, name='sp_res_' + _nm, what='sp resistive probe ' + _nm, set=dict(_s)))
    CASES.append(dict(_CSR, name='cs_res_' + _nm, what='cs resistive probe ' + _nm,
                      set=dict(_CSRS, **_s)))

# ---------------------------------------------------------------- merged lateral cells
# RSG E production: merge_lat_f = 8:4:2 with 2x2 blocks per panel, nghost 3.  A "narrow"
# group (f ceil(2/f) > ng-1, i.e. f >= 4 at ng 3) takes a zero slope at EVERY MeshBlock
# edge (hydro_fluxes.cpp), so the panel midlines become first-order in 2x2 layouts.
_ML = {'mesh/nx1': 8, 'mesh/nx2': 32, 'mesh/nx3': 32, 'mesh/nghost': 3,
       'mesh/cs_seam_blockend': 'extend', 'hydro/merge_lat_rc': 1.0e3}
CASES += [
    var(CS_HYD, 'cs_mlat_f4', 'cs hydro blast, merge_lat_f 4 everywhere (ng 3), extend',
        extra=dict(_ML, **{'hydro/merge_lat_f': 4}), expect='FAIL',
        layouts=[L('b1x1_np1', (8, 32, 32)), L('b2x2_np1', (8, 16, 16)),
                 L('b2x2_np4', (8, 16, 16), 4)]),
    var(CS_HYD, 'cs_mlat_f2', 'cs hydro blast, merge_lat_f 2 (not narrow at ng 3), extend',
        extra=dict(_ML, **{'hydro/merge_lat_f': 2}),
        layouts=[L('b1x1_np1', (8, 32, 32)), L('b2x2_np1', (8, 16, 16)),
                 L('b2x2_np4', (8, 16, 16), 4)]),
]
# radial-BC probes of the resistive failure (located at i = ie on internal block ends)
CASES += [
    dict(_SPR, name='sp_res_ohm_outflow', what='sp Ohmic, OUTFLOW radial BCs',
         set={'mesh/ix1_bc': 'outflow', 'mesh/ox1_bc': 'outflow', 'problem/user_bcs': 'false'}),
    dict(_CSR, name='cs_res_ohm_reflect', what='cs Ohmic, REFLECT radial BCs',
         set=dict(_CSRS, **{'mesh/ix1_bc': 'reflect', 'mesh/ox1_bc': 'reflect',
                            'problem/user_bcs': 'false'})),
]
CASES += [
    dict(_SPR, name='sp_res_ohm_polavg', what='sp Ohmic + use_polar_average_eresist',
         set={'mesh/use_polar_average_eresist': 'true'}),
]

# ---------------------------------------------------------------- red_giant IC seed
# The default vpert seed is periodic over each MeshBlock (O(1) layout dependence of the
# IC); <problem>/vpert_global = true takes it over the mesh.  ncyc 1: the seed itself.
_RG = dict(bin='rg', deck='rg_cs.athinput', ncyc=1, vars=['hydro_u'], timeout=3600,
           layouts=[L('b1x1_np1', (48, 16, 16)), L('b2x2_np4', (48, 8, 8), 4)])
CASES += [
    dict(_RG, name='rg_vpert_block', what='red_giant cs, default vpert seed', expect='FAIL',
         set={'mesh/nx1': 48, 'mesh/cs_seam_blockend': 'extend', 'problem/rt_ck': 'false'}),
    dict(_RG, name='rg_vpert_global', what='red_giant cs, vpert_global = true',
         set={'mesh/nx1': 48, 'mesh/cs_seam_blockend': 'extend', 'problem/rt_ck': 'false',
              'problem/vpert_global': 'true'}),
]

# ---------------------------------------------------------------- exact-solution tests
# layoutaudit-1103: the Ivan et al. (JCP 2013) convergence problems (cs_test iprob 17 hydro
# supersonic expansion, 18 MHD manufactured solution with a user source), C16 x 8, extend.
# Smooth, steady, non-zero flow through every seam and the radial user BCs (exact inflow,
# linearly extrapolated outflow).  Known: MHD + extend fails at > 2 blocks per panel side
# (being fixed elsewhere), so cs_ivan_mhd may FAIL only in b4x4_np6.
_IV = {'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16, 'mesh/cs_seam_blockend': 'extend'}
CASES += [
    dict(name='cs_ivan_hyd', bin='cs', deck='cs_ivan_hyd.athinput', ncyc=10,
         vars=['hydro_u'], what='cs Ivan 2013 supersonic expansion (iprob 17), hllc, extend',
         set=dict(_IV), layouts=cs_layouts(16, 8)),
    dict(name='cs_ivan_mhd', bin='cs', deck='cs_ivan_mhd.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'],
         what='cs Ivan 2013 MHD manufactured solution (iprob 18, user source), hlld, extend',
         set=dict(_IV), layouts=cs_layouts(16, 8)),
]
# separate the rank count from the phi layout (vet_gd off)
_accn = [c for c in CASES if c['name'] == 'acc_sp_novetgd'][0]
CASES.append(dict(_accn, name='acc_sp_novetgd_sep', what='accretor, vet_gd off: rank count '
                  'only (ph2_np4) vs phi layout only (ph4_np2)',
                  layouts=[L('ph2_np2', (500, 48, 32), 2), L('ph2_np4', (500, 48, 32), 4),
                           L('ph4_np2', (500, 48, 16), 2)]))
