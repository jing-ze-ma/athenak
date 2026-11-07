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
        'cs_seam_blockend = clamp (the 11-03 bug, kept visible)', expect='FAIL',
        extra={'mesh/cs_seam_blockend': 'clamp'}),   # explicit: extend is the fresh default since prod1103
    var(CS_HYD, 'cs_hyd_extend', 'cs hydro blast, cs_seam_blockend = extend',
        extra={'mesh/cs_seam_blockend': 'extend'}),
    var(CS_HYD, 'cs_hyd_ng3', 'cs hydro blast, nghost 3 ppm (HG/RSG halo depth), extend',
        extra={'mesh/cs_seam_blockend': 'extend', 'mesh/nghost': 3,
               'hydro/reconstruct': 'ppm4'}),
    var(CS_MHD, 'cs_mhd_clamp', 'cs MHD blast (hlld, CT), clamp', expect='FAIL',
        extra={'mesh/cs_seam_blockend': 'clamp'}),
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
              'mhd/eta_ad': 0.05, 'mhd/use_rkg_sts': 'true',
              'mesh/cs_seam_blockend': 'clamp'},
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
         perturb='problem/b0',
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
    dvar(DHJ_M, 'dhj_r1m_clamp', 'r1m v6 deck, clamp (the old default; explicit since the prod1103 extend default)',
         {'mesh/cs_seam_blockend': 'clamp'}, expect='FAIL'),
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
# the IC field differs between layouts in the last bit (vector potential on per-block
# coordinates): the 1-ulp reference perturbs the field amplitude b0, not the cfl
_SPR = dict(bin='sp', deck='spherical_polar_resist.athinput', ncyc=10, vars=['mhd_u', 'mhd_bcc'],
            perturb='problem/b0',
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
    var(CS_HYD, 'cs_mlat_f4', 'cs hydro blast, merge_lat_f 4 everywhere (ng 3), extend, '
        'merge_lat_blockend = flat (the old rule, kept visible)',
        extra=dict(_ML, **{'hydro/merge_lat_f': 4, 'hydro/merge_lat_blockend': 'flat'}),
        expect='FAIL',
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

# merge_lat_blockend = exchange (layoutaudit-1103): the gates
_MLL = [L('b1x1_np1', (8, 32, 32)), L('b2x2_np1', (8, 16, 16)), L('b2x2_np4', (8, 16, 16), 4),
        L('b2x2_np6', (8, 16, 16), 6), L('b1x1_np6', (8, 32, 32), 6)]
CASES += [
    var(CS_HYD, 'cs_mlat_f4x', 'cs hydro blast, merge_lat_f 4, merge_lat_blockend = exchange',
        extra=dict(_ML, **{'hydro/merge_lat_f': 4, 'hydro/merge_lat_blockend': 'exchange'}),
        layouts=_MLL),
    var(CS_HYD, 'cs_mlat_f8x', 'cs hydro blast, merge_lat_f 8, merge_lat_blockend = exchange',
        extra=dict(_ML, **{'hydro/merge_lat_f': 8, 'hydro/merge_lat_blockend': 'exchange'}),
        layouts=_MLL),
    var(CS_HYD, 'cs_mlat_842x', 'cs hydro blast, merge_lat_f 8:4:2 shells (RSG pattern), '
        'exchange', extra=dict(_ML, **{'hydro/merge_lat_f': '8:4:2',
                                       'hydro/merge_lat_rc': '1.3:1.6:1.85',
                                       'hydro/merge_lat_blockend': 'exchange'}),
        layouts=_MLL),
    var(CS_HYD, 'cs_mlat_f2x', 'cs hydro blast, merge_lat_f 2 (not narrow), exchange',
        extra=dict(_ML, **{'hydro/merge_lat_f': 2, 'hydro/merge_lat_blockend': 'exchange'}),
        layouts=_MLL),
]
CASES += [
    var(CS_HYD, 'cs_mlat_rsgkeys', 'cs hydro blast, RSG E key set: merge_lat 8:4:2, ng 3, '
        'cs_seam_flux shared, fofc, exchange',
        extra=dict(_ML, **{'hydro/merge_lat_f': '8:4:2', 'hydro/merge_lat_rc': '1.3:1.6:1.85',
                           'hydro/merge_lat_blockend': 'exchange', 'hydro/fofc': 'true',
                           'mesh/cs_seam_flux': 'shared'}),
        layouts=_MLL),
    var(CS_HYD, 'cs_mlat_rsgkeys_flat', 'as cs_mlat_rsgkeys with merge_lat_blockend = flat '
        '(= RSG E production now)',
        extra=dict(_ML, **{'hydro/merge_lat_f': '8:4:2', 'hydro/merge_lat_rc': '1.3:1.6:1.85',
                           'hydro/merge_lat_blockend': 'flat', 'hydro/fofc': 'true',
                           'mesh/cs_seam_flux': 'shared'}),
        expect='FAIL', layouts=_MLL[:3]),
]

# RSG E production deck (he_star_m1, merge_lat 8:4:2, ng 3, cs_seam_flux shared, fofc, M1),
# laterally C32 (2x2 blocks of 16 = 2 x f 8), the production radial grid (nx1 548).
# perturb he_gm: dt may be set by the vet_gd rotation alignment, not the cfl.
_RSG = dict(bin='he', deck='he_rsg.athinput', ncyc=3, timeout=14400, perturb='problem/he_gm',
            vars=['hydro_u'],
            layouts=[L('b1x1_np6', (548, 32, 32), 6), L('b2x2_np6', (548, 16, 16), 6)])
_RSGS = {'mesh/nx2': 32, 'mesh/nx3': 32, 'mesh/cs_seam_blockend': 'extend',
         'hydro/merge_lat_blockend': 'exchange'}
# share-1104 (11-05): split.  rsg_mlat_x = the PRODUCTION check (merge_lat_blockend =
# exchange, RSG E since 11-03); rsg_mlat_flat = the old flat rule, an expected FAIL (narrow
# merge groups take a zero slope at internal block ends: 1e-4..1e-3 in mom at i 1-5).
# Known on 6a11ee57e (share_1104/RSG_MLAT_LAYOUT.md): the exchange layouts differ ~100x the
# 1-ulp reference in the front band i 317-323 only, seeded by a layout-dependent IC (whole
# shells at c0) plus a one-time jump at the block-end corners in cycle 1, flat afterwards.
CASES += [
    dict(_RSG, name='rsg_mlat_x', what='RSG E deck C32, merge_lat_blockend = exchange '
         '(the production rule)', set=dict(_RSGS)),
    dict(_RSG, name='rsg_mlat_flat', what='RSG E deck C32, merge_lat_blockend = flat (old rule)',
         expect='FAIL', set=dict(_RSGS),
         layouts=[L('b1x1_np6', (548, 32, 32), 6),
                  L('b2x2flat_np6', (548, 16, 16), 6, set={'hydro/merge_lat_blockend': 'flat'})]),
]
CASES += [
    dict(_SPR, name='sp_res_ohm_nopolavg', what='sp Ohmic, use_polar_average_eresist = false',
         ncyc=1, set={'mesh/use_polar_average_eresist': 'false'}),
    dict(_SPR, name='sp_res_ohm_1c', what='sp Ohmic, 1 cycle', ncyc=1, set={}),
]

# ---------------------------------------------------------------- cube-vertex corner probe
# (USER 11-03) does the EXACT r1m v8 configuration read the cell-centred cube-vertex
# corner block?  Binary family dhjvtx = deep_hot_jupiter_rt built from vertexacc-1103
# (<mesh>/cs_vertex_cc_perturb = eps multiplies the extrapolated corner block by 1+eps).
# "Layouts" here are the same layout with the perturbation: BITWISE = never read.
_V8 = dict(bin='dhjvtx', deck='dhj_r1m_v8.athinput', ncyc=10, vars=['mhd_u', 'mhd_bcc'],
           links=DHJ_LINKS, timeout=7200,
           set={'mesh/nx1': 32, 'mesh/nx2': 8, 'mesh/nx3': 8, 'time/ndiag': 1},
           layouts=[L('ref', (32, 8, 8)), L('vp1e-2', (32, 8, 8),
                                            set={'mesh/cs_vertex_cc_perturb': 1e-2}),
                    L('vp1e200', (32, 8, 8), set={'mesh/cs_vertex_cc_perturb': 1e200}),
                    L('vp1e-2_np6', (32, 8, 8), 6, set={'mesh/cs_vertex_cc_perturb': 1e-2})])
CASES += [dict(_V8, name='r1m_v8_vertex_cc', what='r1m v8 deck C8x32 (1 block/panel), corner '
               'block x (1+eps): BITWISE = the corner block is never read')]
# sp pole-row probes (the sp_resist failure is already in the IDEAL corner EMF e2 at j = js)
_SPI = dict(_SPR, ncyc=1)
CASES += [
    dict(_SPI, name='sp_ideal_res_deck', what='sp_resist deck, NO resistivity',
         set={'mhd/eta_ohm_const': 0.0}),
    dict(_SPI, name='sp_ideal_nodiss', what='no resistivity, polar_emf_diss = false',
         set={'mhd/eta_ohm_const': 0.0, 'mhd/polar_emf_diss': 'false'}),
    dict(_SPI, name='sp_ideal_gs07', what='no resistivity, emf_scheme = gs07',
         set={'mhd/eta_ohm_const': 0.0, 'mhd/emf_scheme': 'gs07'}),
    dict(_SPI, name='sp_ideal_gs07_nodiss', what='no resistivity, gs07, no polar_emf_diss',
         set={'mhd/eta_ohm_const': 0.0, 'mhd/emf_scheme': 'gs07', 'mhd/polar_emf_diss': 'false'}),
]
# the IC field differs between layouts in the last bit (vector potential on per-block
# coordinates), so the reference for these is a 1-ulp IC perturbation of the field
CASES += [
    dict(_SPI, name='sp_ideal_gs07_iculp', what='no resistivity, gs07; 1-ulp IC (b0) reference',
         perturb='problem/b0', set={'mhd/eta_ohm_const': 0.0, 'mhd/emf_scheme': 'gs07'}),
    dict(_SPI, name='sp_ideal_uct_iculp', what='no resistivity, uct_hlld + polar_emf_diss; '
         '1-ulp IC (b0) reference', perturb='problem/b0', set={'mhd/eta_ohm_const': 0.0}),
    dict(_SPI, name='sp_res_ohm_iculp', what='Ohmic; 1-ulp IC (b0) reference',
         perturb='problem/b0', set={}),
]

# ---------------------------------------------------------------- fcorder3-1103
# <mesh>/cs_seam_fc_order = 3 (third-order seam ghost face fields): must stay ROUNDOFF under
# cs_seam_blockend = extend at every layout (its along stencils are clamped to within ng of
# the sample, so no layout reads beyond its extended bounds).
_F3 = {'mesh/cs_seam_fc_order': 3}
CASES += [
    var(CS_MHD, 'cs_mhd_ext_fc3', 'cs MHD blast, extend, cs_seam_fc_order 3',
        extra=dict(EXT, **_F3)),
    var(CS_MHD, 'cs_mhd_ext_n32_fc3', 'cs MHD blast C32 1 vs 4x4 / 2x2 / 8x8, extend, fc order 3',
        extra=dict(EXT, **_F3, **{'mesh/nx2': 32, 'mesh/nx3': 32}),
        layouts=[L('b1x1_np1', (8, 32, 32)), L('b4x4_np1', (8, 8, 8)),
                 L('b2x2_np1', (8, 16, 16)), L('b8x8_np6', (8, 4, 4), 6)]),
    var(CS_MHD, 'cs_mhd_ng3_fc3', 'cs MHD blast, nghost 3, extend, fc order 3',
        extra=dict(EXT, **_F3, **{'mesh/nghost': 3})),
    dict(name='cs_resist_rkg_fc3', bin='cs', deck='cubed_sphere_resist.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'], what='cs Ohmic + ambipolar, RKG STS, extend, fc order 3',
         set={'mesh/nx1': 8, 'mesh/nx2': 16, 'mesh/nx3': 16,
              'mesh/cs_seam_blockend': 'extend', 'mhd/ambipolar': 'constant',
              'mhd/eta_ad': 0.05, 'mhd/use_rkg_sts': 'true', 'mesh/cs_seam_fc_order': 3},
         layouts=cs_layouts(16, 8)),
    dict(name='cs_ivan_mhd_fc3', bin='cs', deck='cs_ivan_mhd.athinput', ncyc=10,
         vars=['mhd_u', 'mhd_bcc'], what='cs Ivan MHD (iprob 18), extend, fc order 3',
         set=dict(_IV, **_F3), layouts=cs_layouts(16, 8)),
    dvar(DHJ_M, 'dhj_r1m_fc3', 'r1m v6 MHD deck, extend, cs_seam_fc_order 3', dict(EXT, **_F3)),
    dict(_V8, name='r1m_v8_vertex_cc_fc3', what='r1m v8 deck C8x32, corner block x (1+eps), '
         'fc order 3: BITWISE = the corner block is never read', set=dict(_V8['set'], **_F3)),
]
for _nm, _s in (('ohm', {}), ('ohm_rkg', {'mhd/use_rkg_sts': 'true'}),
                ('amb', {'mhd/ambipolar': 'constant', 'mhd/eta_ad': 0.05,
                         'mhd/ohmic_resistivity': 'none'}),
                ('amb_rkg', {'mhd/ambipolar': 'constant', 'mhd/eta_ad': 0.05,
                             'mhd/ohmic_resistivity': 'none', 'mhd/use_rkg_sts': 'true'})):
    CASES.append(dict(_CSR, name='cs_res_' + _nm + '_fc3', what='cs resistive probe ' + _nm
                      + ', fc order 3', set=dict(_CSRS, **_s, **_F3)))

# ---------------------------------------------------------------- fc3blk-1103
# order 3 with 4-cell blocks (r1m v6 deck C8x32, 1x1 vs 2x2, cs_vertex_3cell off).  The
# cfl 1-ulp reference reaches the stiff ck operator only weakly; the second case takes a
# 1-ulp IC field perturbation (problem/bbot_gauss) as the reference, which is the honest
# round-off spread of this deck after one ck step (1e-11 in ener at either order).
_R4L = [L('b1x1_np1', (32, 8, 8)), L('b2x2_np1', (32, 4, 4))]
_R4S = dict(EXT, **{'mesh/cs_vertex_3cell': 'false', 'mesh/cs_seam_fc_order': 3,
                    'mesh/nx2': 8, 'mesh/nx3': 8})
CASES += [
    dict(dvar(DHJ_M, 'r1m4_fc3', 'r1m v6 C8 4-cell blocks, order 3, cfl 1-ulp reference', _R4S),
         layouts=_R4L),
    dict(dvar(DHJ_M, 'r1m4_fc3_iculp', 'r1m v6 C8 4-cell blocks, order 3, IC field 1-ulp '
              'reference', _R4S), layouts=_R4L, perturb='problem/bbot_gauss'),
    dict(dvar(DHJ_M, 'r1m4_fc2_iculp', 'as r1m4_fc3_iculp, order 2',
              dict(_R4S, **{'mesh/cs_seam_fc_order': 2})), layouts=_R4L,
         perturb='problem/bbot_gauss'),
]

# ---------------------------------------------------------------- layinv-1105
# <mesh>/layout_invariant (fresh default): the accretor full sphere with a THETA split of
# the stretched (plateau) theta grid as well as phi splits, ranks 2..4.  The accretor's
# M1 preconditioner is mg (block-local): bitwise is NOT expected after cycle 0 (USER
# 11-05 option a: mg/rbgs decks are gated at the tolerance level, see README).
CASES.append(dict(_acc, name='acc_sp_th', what='accretor deck, theta + phi splits '
                  '(stretched plateau theta), layout_invariant; (a1) gate: mg is block-local, '
                  'PASS = no worse than the old binary (family acc_old) on the same layouts, '
                  'with acc_sp_th_line bitwise', ref_old=True,
                  layouts=[L('ph2_np2', (500, 48, 32), 2), L('th2ph2_np4', (500, 24, 32), 4),
                           L('ph4_np4', (500, 48, 16), 4)]))
CASES.append(dict(_acc, name='acc_sp_th_line', what='accretor deck, theta + phi splits, '
                  'implicit_precond = line (layout-independent preconditioner: bitwise '
                  'expected)', set=dict(_acc['set'], **{'rad_m1/implicit_precond': 'line'}),
                  require_bitwise=True,
                  layouts=[L('ph2_np2', (500, 48, 32), 2), L('th2ph2_np4', (500, 24, 32), 4),
                           L('ph4_np4', (500, 48, 16), 4)]))

# layinv-1105: the BSG column deck (bsgLI_exact derivative: sp WEDGE pi/3 x pi/3, 256 x 8 x 8,
# vet_gd twin trans + symref col, implicit_precond = line, Linhao IC with a 5 % lateral
# perturbation; needs the linhao-prod-1104 IC, i.e. prod-1102 after the athpp-doface-1105
# merge).  layout_invariant on (fresh default): BITWISE required across theta, phi and
# theta+phi splits and ranks.
_HCG = dict(bin='he', deck='he_col_gd.athinput', ncyc=3, vars=['hydro_w', 'm1'],
            set={'mesh/nx2': 8, 'mesh/nx3': 8, 'problem/he_ic_linhao_amp': 0.05},
            require_bitwise=True)
CASES += [dict(_HCG, name='he_col_gd', what='BSG col_gd_df deck 256 x 8 x 8, layout_invariant',
               layouts=[L('b1_np1', (256, 8, 8)), L('ph2_np1', (256, 8, 4)),
                        L('th2_np1', (256, 4, 8)), L('th2ph2_np1', (256, 4, 4)),
                        L('th2ph2_np4', (256, 4, 4), 4)]),
          dict(_HCG, name='he_col_gd_off', what='the same deck with layout_invariant = false '
               '(the old block-ordered paths: expected to differ)', expect='FAIL',
               require_bitwise=True,
               set=dict(_HCG['set'], **{'mesh/layout_invariant': 'false'}),
               layouts=[L('b1_np1', (256, 8, 8)), L('th2ph2_np1', (256, 4, 4))])]

# layinv-1105 (a1, USER 11-05): the M1 decks with a block-local preconditioner (mg) are gated
# against the OLD binary's same layouts (families he_old / acc_old in the bins file)
for _c in CASES:
    if _c['name'] in ('he_hg', 'rsg_mlat_x', 'acc_sp', 'acc_sp_novetgd', 'acc_sp_novetgd_sep'):
        _c['ref_old'] = True

# wallwb-1105: HG deck (closed inner wall + x1 WB, cubed sphere) with <hydro>/wall_wb_onesided
# = true.  The key touches only the wall cell of each radial column (block-local, no
# neighbour data): layout pairs must be no worse than the OLD binary's (which ignores the
# key; mg preconditioner = block-local, hence ref_old as he_hg).
_hg = [c for c in CASES if c['name'] == 'he_hg'][0]
CASES.append(dict(_hg, name='he_hg_wallwb', ref_old=True,
                  what='HG deck with hydro/wall_wb_onesided = true',
                  set=dict(_hg['set'], **{'hydro/wall_wb_onesided': 'true'}),
                  layouts=[L('b1x1_np3', (704, 16, 16), 3), L('b2x2_np4', (704, 8, 8), 4)]))

# wallwb-1105 v2: the R136a1 column (sp wedge, closed wall, x1 WB, general EOS, M1 vet_col +
# vet_gd), 831 x 8 x 8, seed 0, implicit_precond = line (layout-independent), the fresh-start
# default wall_wb_onesided = true active: theta + phi split and 4 ranks.  ref_old as he_hg.
CASES.append(dict(name='he_r136col', bin='he', deck='he_r136col.athinput', ncyc=3, timeout=3000,
                  perturb='problem/he_gm', vars=['hydro_u', 'm1'], ref_old=True,
                  what='R136a1 column 831 x 8 x 8, line precond, default wall_wb_onesided',
                  set={'rad_m1/implicit_precond': 'line'},
                  layouts=[L('b1_np1', (831, 8, 8)), L('th2ph2_np4', (831, 4, 4), 4)]))

# layoutfix-1106: separate the two contributions to the he_hg layout pair.
# he_hg_noic*: identical analytic column IC (no ic_rad, no he_ic_equil): cycle 0 must be bitwise and
# the pair measures the time step only.  he_hg_line*: implicit_precond = line (layout-independent
# preconditioner): the deck must be bitwise, with and without the IC relaxation.
_hg4 = [L('b1x1_np4', (704, 16, 16), 4), L('b2x2_np4', (704, 8, 8), 4)]
_noic = {'rad_m1/ic_rad': 'none', 'problem/he_ic_equil': 'false'}
_line = {'rad_m1/implicit_precond': 'line', 'rad_m1/cs_precond_line': 'true'}
CASES.append(dict(_hg, name='he_hg_noic', ref_old=False, layouts=_hg4,
                  what=_hg['what'] + ', no IC relaxation', set=dict(_hg['set'], **_noic)))
CASES.append(dict(_hg, name='he_hg_noic_line', ref_old=False, require_bitwise=True, layouts=_hg4,
                  what=_hg['what'] + ', no IC relaxation, line precond',
                  set=dict(_hg['set'], **_noic, **_line)))
CASES.append(dict(_hg, name='he_hg_line', ref_old=False, require_bitwise=True, layouts=_hg4,
                  what=_hg['what'] + ', line precond',
                  set=dict(_hg['set'], **_line)))
# tolerances x 1e-3: does the mg pair fall with the solver tolerance?
_t3 = {'rad_m1/implicit_lin_tol': '1.0e-13', 'rad_m1/implicit_lin_cnorm': '1.0e-13',
       'rad_m1/implicit_tol': '1.0e-11'}
CASES.append(dict(_hg, name='he_hg_noic_t3', ref_old=False, layouts=_hg4,
                  what=_hg['what'] + ', no IC relaxation, tolerances x 1e-3',
                  set=dict(_hg['set'], **_noic, **_t3)))
# Eisenstat-Walker and the relative Picard stop off (the deck's real stops are ew_max 1e-2 and
# tol_rel 0.1, capped at 1e-4): t4 with the deck's tolerances, t5 with them x 1e-3
_ewoff = {'rad_m1/implicit_lin_ew_max': '0.0', 'rad_m1/implicit_tol_rel': '0.0'}
CASES.append(dict(_hg, name='he_hg_noic_t4', ref_old=False, layouts=_hg4,
                  what=_hg['what'] + ', no IC relaxation, EW and tol_rel off',
                  set=dict(_hg['set'], **_noic, **_ewoff)))
CASES.append(dict(_hg, name='he_hg_noic_t5', ref_old=False, layouts=_hg4,
                  what=_hg['what'] + ', no IC relaxation, EW and tol_rel off, tolerances x 1e-3',
                  set=dict(_hg['set'], **_noic, **_ewoff, **_t3)))
# layoutfix-1106: <mesh>/cs_seam_pos_global: the line arm must be BITWISE across layouts
_pg = {'mesh/cs_seam_pos_global': 'true'}
CASES.append(dict(_hg, name='he_hg_noic_line_pg', ref_old=False, require_bitwise=True, layouts=_hg4,
                  what=_hg['what'] + ', no IC relaxation, line precond, cs_seam_pos_global',
                  set=dict(_hg['set'], **_noic, **_line, **_pg)))
CASES.append(dict(_hg, name='he_hg_line_pg', ref_old=False, require_bitwise=True, layouts=_hg4,
                  what=_hg['what'] + ', line precond incl. the IC relaxation, cs_seam_pos_global',
                  set=dict(_hg['set'], **_line, **_pg)))
# the production mg deck with a layout-independent initial state: cycle 0 must be bitwise
_icl = {'rad_m1/ic_rad_precond': 'line'}
CASES.append(dict(_hg, name='he_hg_icline_pg', ref_old=False, layouts=_hg4,
                  what=_hg['what'] + ', mg, IC relaxation with line, cs_seam_pos_global',
                  set=dict(_hg['set'], **_icl, **_pg)))
