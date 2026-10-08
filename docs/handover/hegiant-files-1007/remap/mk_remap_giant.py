#!/usr/bin/env python3
"""he_remap_giant.py = he_mltpp_1002/w256/he_remap_rst.py with the radial potential taken EXACTLY:
old cells Phi = (E - KE - e_int)/rho from the file (per cell), new cells Phi = phi_code.phi_cells (the he_gm_column
potential the pgen rebuilds on the new grid, validated to 1e-14 on the scout rst).  The source assumed a point
mass (Phi = a + b/r fit), wrong for he_gm_column."""
src = open('/viper/ptmp2/jinma/he_giant_1006/remap/he_remap_rst_w256.py').read()
s = src


def rep(o, n):
    global s
    assert s.count(o) == 1, o[:60]
    s = s.replace(o, n)


rep('"""he_remap_rst.py', '"""he_remap_giant.py (he_giant_1006 copy of he_mltpp_1002/w256/he_remap_rst.py, made by '
    'mk_remap_giant.py: the radial path uses the EXACT he_gm_column potential, old from the file, new from '
    'phi_code.py)\nhe_remap_rst.py')
rep("""    phi = (u[4] - ke - G['eint'][0]) / rho
    pm = np.median(phi, axis=(0, 1))
    Am = np.vstack([np.ones_like(xo), 1.0 / xo]).T
    cf = np.linalg.lstsq(Am[ao], pm[ao], rcond=None)[0]
    st['phi_fit'] = (cf, float(np.abs((Am @ cf - pm) / pm)[ao].max()))
    phin = cf[0] + cf[1] / xn
""", """    sys.path.insert(0, '/viper/ptmp2/jinma/he_giant_1006/remap')
    from phi_code import phi_cells
    gm_ = float(get_param(r.text, 'problem', 'he_gm'))
    icf_ = get_param(r.text, 'problem', 'he_ic_file')
    nf_ = int(get_param(r.text, 'problem', 'he_nfine'))
    _, pho = phi_cells(fo, ng, icf_, gm_, nf_)
    _, phin = phi_cells(fn, ng, icf_, gm_, nf_)
    if 'eint' in G:
        phi = (u[4] - ke - G['eint'][0]) / rho
        st['phi_old_file_vs_code'] = float((np.abs(phi - pho)[..., ao]).max() / np.abs(pho).max())
    else:       # no e_int cache (a file written by this script): the code's Phi (1e-14) and e_int = U - KE
        phi = np.broadcast_to(pho, rho.shape)
        G['eint'] = (u[4] - ke - rho * phi)[None]
        st['phi_old'] = 'code (no e_int cache in the file)'
""")
rep("for nm, qo, qn in (('mass', rho, hyd[0]), ('U', U, Un), ('m_r', u[1], hyd[1])):",
    "for nm, qo, qn in (('mass', rho, hyd[0]), ('U', U, Un), ('m_r', u[1], hyd[1]),\n"
    "                       ('r_m_th', u[2] * xo, hyd[2] * xn), ('r_m_ph', u[3] * xo, hyd[3] * xn)):")
rep("""    if keep:
        keepset |= {'wt', 'wd', 'pred', 'eint', 'rss'}""", """    if keep:
        keepset |= {'wt', 'wd', 'pred', 'eint', 'rss'}
    elif KEEP_M1:             # --keep-m1-caches: the M1 step predictor and rad_signal_speed inputs only
        keepset |= {'pred', 'rss'}""")
rep("""        if nm in dropped and not keep:
            continue""", """        if nm in dropped and not keep and not (KEEP_M1 and nm in ('pred', 'rss')):
            continue""")
rep("""    ap.add_argument('--keep-dt', action='store_true')""", """    ap.add_argument('--keep-dt', action='store_true')
    ap.add_argument('--keep-m1-caches', action='store_true',
                    help='keep only the M1 predictor + rad_signal_speed blocks (EOS caches dropped)')
    ap.add_argument('--m1-conservative', action='store_true',
                    help='old radial M1 remap: E and F conservative (breaks free streaming)')""")
rep("""    a = ap.parse_args()""", """    a = ap.parse_args()
    global KEEP_M1, M1_RATIO
    KEEP_M1 = a.keep_m1_caches
    M1_RATIO = not a.m1_conservative""")
rep("""NIND = 19""", """KEEP_M1 = False
M1_RATIO = True   # he_giant 10-07: radiation F remapped as the flux RATIO (see radial_remap); --m1-conservative
NIND = 19""")
rep("""    qs = [w[0], w[1], w[2] * xo, w[3] * xo]
    ss = [slope_of(q) for q in qs]
    zero = np.zeros(w[0][..., ao].shape, dtype=bool)
    for it in range(6):""", """    if M1_RATIO:
        # he_giant fix (10-07): E conservative (limited linear, slope dropped where E <= 0); F = c E_new x the
        # E-weighted overlap mean of the old flux RATIO f = F/(cE), with |f| set to the E-weighted mean of |f_old|.
        # Free streaming (|f| = 1 to round-off in the thin region) survives; a conservative F remap broke it
        # (|F|/cE 0.93-0.9997 above 65 Rsun) and the implicit M1 then needed ~110 Picard passes per solve.
        zero = np.zeros(w[0][..., ao].shape, dtype=bool)
        se = slope_of(w[0])
        for it in range(6):
            e0 = remap1(w[0], se, zero)
            badn = (e0 <= 0.0)[..., an]
            if not badn.any():
                break
            zero |= old_of(badn)
        st['rad_m1_E_slope_dropped_old_cells'] = int(zero.sum())
        Ea = w[0][..., ao]
        rat = [w[n][..., ao] / (cl * Ea) for n in (1, 2, 3)]
        mag = np.sqrt(rat[0]**2 + rat[1]**2 + rat[2]**2)
        wsum = Ea @ A.T
        avg = [(q * Ea) @ A.T / wsum for q in rat]
        mavg = np.minimum((mag * Ea) @ A.T / wsum, 1.0)
        aabs = np.sqrt(avg[0]**2 + avg[1]**2 + avg[2]**2)
        sc = np.where(aabs > 1e-300, mavg / np.maximum(aabs, 1e-300), 0.0)
        m1 = np.empty((4,) + e0.shape)
        m1[0] = e0
        for n in (1, 2, 3):
            m1[n][..., an] = cl * e0[..., an] * avg[n - 1] * sc
            m1[n][..., :ng] = w[n][..., :ng]
            m1[n][..., ng + n1new:] = w[n][..., ng + n1:]
        to, tn = tot(w[0], dVo, ao), tot(m1[0], dVn, an)
        st['rad_cons_Erad'] = float((np.abs(tn - to) / to).max())
        to, tn = tot(w[1], dVo, ao), tot(m1[1], dVn, an)
        st['rad_Fr_change_rel'] = float((np.abs(tn - to) / np.abs(to).max()).max())
        st['rad_F_ratio_max'] = float((np.sqrt(m1[1]**2 + m1[2]**2 + m1[3]**2) / (cl * m1[0]))[..., an].max())
        Gn['m1'] = m1
    qs = [w[0], w[1], w[2] * xo, w[3] * xo]
    ss = [slope_of(q) for q in qs] if not M1_RATIO else []
    zero = np.zeros(w[0][..., ao].shape, dtype=bool)
    for it in range(6 if not M1_RATIO else 0):""")
rep("""    m1 = np.stack([e0, f1, f2, f3])
    to, tn = tot(w[0], dVo, ao), tot(m1[0], dVn, an)""", """    if M1_RATIO:
        pass
    else:
      m1 = np.stack([e0, f1, f2, f3])
    to, tn = tot(w[0], dVo, ao), tot(m1[0], dVn, an)""")
rep("""    near = np.abs(xn[:, None] - xo[None, :]).argmin(axis=1)""", """    if ZERO_FACES:
        # he_giant fix 2 (10-07): the stored M1 face fluxes f0x1..3 seed the lagged closure ratio r0 (each face
        # normalised by its UPWIND cell E, rad_m1_impl_lag_kernel.hpp); interpolated onto new faces they are not
        # consistent with the remapped cells (in the thin region f1/(cE) was off by O(1)).  Zero = the code's own
        # cold-start path: r0 from the cell flux on the first pass, faces recomputed at the end of every pass.
        for nm in ('f1', 'f2', 'f3'):
            Gn[nm] = np.zeros_like(Gn[nm])
        st['m1_faces'] = 'zeroed (cold-start closure lag from the cell F)'
    near = np.abs(xn[:, None] - xo[None, :]).argmin(axis=1)""")
rep("""M1_RATIO = True""", """M1_RATIO = True
ZERO_FACES = True  # he_giant 10-07 fix 2: zero the stored M1 face fluxes on a radial remap; --keep-faces""")
rep("""    global KEEP_M1, M1_RATIO""", """    global KEEP_M1, M1_RATIO, ZERO_FACES
    ZERO_FACES = not a.keep_faces""")
rep("""    ap.add_argument('--m1-conservative', action='store_true',""", """    ap.add_argument('--keep-faces', action='store_true',
                    help='radial remap: interpolate the M1 face fluxes instead of zeroing them')
    ap.add_argument('--m1-conservative', action='store_true',""")
open('/viper/ptmp2/jinma/he_giant_1006/remap/he_remap_giant.py', 'w').write(s)
print('ok')
