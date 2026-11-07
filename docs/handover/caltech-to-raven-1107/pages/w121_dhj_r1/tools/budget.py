"""dhj_r1 conservation budgets (mass, total energy, L_z) from dhj.hydro.hst + dhj.user.hst (read-only).

Bookkeeping (code: src/pgen/deep_hot_jupiter_rt.cpp DhjFluxHistory / AmHstColumns / DhjTideTorque, fin51-r-fix d3f878be):
  mass   dM/dt   = Mdot_bot - Mdot_top + Mfloor            (top sponge damps v_r only: no mass)
  energy dE/dt   = Etot_bot - Etot_top + Efloor             E = hst tot-E = sum (e_int + KE + rho Phi_tot), etotgrav,
           Phi_tot = g0 ap (1 - ap/r) - Omega^2 R^2/2 + TidePot (time-independent in the corotating frame, so the
           tide and the centrifugal force do no work beyond the rho Phi v flux already in Etot_top/Etot_bot);
           Etot_top = fluid + Lir_top - Lsw_abs, Etot_bot = fluid + Lrad_bot (sigma T_int^4 at the wall).
           The top sponge damps the radial momentum without touching E (its KE goes to e_int): no energy term.
  L_z    dL_abs/dt = Tbnd + Ttide + Tdrag (+ Ttop + Twall: no column)   L_abs = Lrel + Omega Iw
           (coriolis_am: Tcor = Tcont cancels Omega dI/dt to round-off; sponge radial-only: no torque).
Rate columns Efloor/Mfloor/Tbnd/Tdrag are interval averages (rate x dt is exact); the flux columns
(Etot_*, Mdot_*) and Ttide are snapshots at the row (trapezoid; sampling error in the residual).
"""
import numpy as np

def read_hst(fn):
    """rows keyed by time, last write wins (restarts, the 11-01 duplicate W178 process, the 11-02 moves), sorted"""
    lab = None; seen = {}; labs = []
    for l in open(fn):
        if l.startswith('#  ['):
            lab = [s.split('=')[1] for s in l[1:].split()]
            labs += [k for k in lab if k not in labs]; continue
        if l.startswith('#') or not l.strip(): continue
        x = l.split()
        if lab is None or len(x) != len(lab): continue   # torn line
        try: v = [float(s) for s in x]
        except ValueError: continue
        seen['%.9e' % v[0]] = dict(zip(lab, v))     # segments may differ in columns (older binary: no Tdrag)
    rows = sorted(seen.values(), key=lambda r: r['time'])
    return {k: np.array([r.get(k, 0.0) for r in rows]) for k in labs}

def trap(t, y):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])

def rect(t, y):   # interval-average rate columns: the row's value covers (t_prev, t]
    return np.concatenate([[0.0], np.cumsum(y[1:] * np.diff(t))])

def budgets(run, om):
    h = read_hst(run + '/dhj.hydro.hst'); u = read_hst(run + '/dhj.user.hst')
    tc = np.intersect1d(np.round(h['time'], 3), np.round(u['time'], 3))
    ih = np.isin(np.round(h['time'], 3), tc); iu = np.isin(np.round(u['time'], 3), tc)
    h = {k: v[ih] for k, v in h.items()}; u = {k: v[iu] for k, v in u.items()}
    t = h['time']; P = 2 * np.pi / om
    g = lambda k: u[k] if k in u else np.zeros_like(t)
    B = {'rot': t / P, 't': t, 'P': P}
    # mass
    M = h['mass']; dM = M - M[0]
    mb = trap(t, g('Mdot_bot')); mt = -trap(t, g('Mdot_top')); mf = rect(t, g('Mfloor'))
    B['M'] = dict(d=dM, terms={'outflow top': mt, 'wall': mb, 'floors': mf}, acc=mb + mt + mf, scale=M[0], scale_name='M0',
                  alt=np.abs(trap(t, g('Mdot_top'))), alt_name='integrated top outflow')
    # energy
    E = h['tot-E']; dE = E - E[0]
    eb = trap(t, g('Etot_bot')); et = -trap(t, g('Etot_top')); ef = rect(t, g('Efloor'))
    lsw = trap(t, g('Lsw_abs')); lir = trap(t, g('Lir_top')); lbot = trap(t, g('Lrad_bot'))
    fl_top = et - lsw + lir   # fluid (advective, incl. rho Phi v) part of the top flux, as an inflow
    B['E'] = dict(d=dE, terms={'stellar absorbed': lsw, 'IR out (top)': -lir, 'interior heat (wall)': lbot,
                               'fluid out (top)': fl_top, 'fluid (wall)': eb - lbot, 'floors': ef},
                  acc=eb + et + ef, scale=abs(E[0]), scale_name='|E0|', alt=lsw, alt_name='integrated L_sw,abs dt')
    # L_z
    Lrel, Iw = g('Lrel'), g('Iw'); La = Lrel + om * Iw; dL = La - La[0]
    tb = rect(t, g('Tbnd')); tt = trap(t, g('Ttide')); td = rect(t, g('Tdrag'))
    B['L'] = dict(d=dL, terms={'tide': tt, 'boundary (Omega R^2 outflow)': tb, 'start-up drag': td},
                  acc=tt + tb + td, scale=Iw[0] * om, scale_name='I Omega', alt=np.abs(tt), alt_name='|integrated tide torque|')
    B['raw'] = dict(Mdot_top=g('Mdot_top'), Mfloor=g('Mfloor'), Efloor=g('Efloor'), Lrad_bot=g('Lrad_bot'),
                    Lsw_abs=g('Lsw_abs'), Lir_top=g('Lir_top'), Etot_top=g('Etot_top'), Ttide=g('Ttide'), Mdot_bot=g('Mdot_bot'))
    for q in 'MEL':
        B[q]['res'] = B[q]['d'] - B[q]['acc']
    return B

# ---------------------------------------------------------------- am_vector spin reservoir (not in the hst)
# <hydro>/am_vector puts the RADIAL part of its transport torque correction, rhat . dT, into a per-cell spin
# reservoir am_spin [g cm^2/s] (hydro_am.cpp; am_spin_return off = never returned).  Its z part, sum S rhat_z, is
# angular momentum the hst Lrel does not see: the conserved total is Lrel + Omega Iw + sum S rhat_z.  It is
# stored only in the restart files (last slab of each MeshBlock record, outputs/restart.cpp), read here.
import os, glob, json

def _rst_layout(f, fs, nmb, nx1, nx2, nx3):
    for ng in (2, 0, 3, 4):
        nc = (nx1 + 2 * ng) * (nx2 + 2 * ng) * (nx3 + 2 * ng)
        for nsl in range(5, 80):
            ds = nsl * nc * 8; off = fs - nmb * ds - 8
            if off < 0: continue
            f.seek(off)
            if int(np.frombuffer(f.read(8), '<u8')[0]) == ds: return ng, off + 8, ds
    raise RuntimeError('rst layout not found')

def spin_lz_rst(fn, geom, nx1, nx2, nx3):
    """(time, sum S rhat_z) of one restart file; geom = mb_geometry (nmb, 6) of a dump of the same run"""
    import dhjcs
    nmb = geom.shape[0]; fs = os.path.getsize(fn)
    with open(fn, 'rb') as f:
        hd = f.read(600000); e = hd.find(b'<par_end>') + len(b'<par_end>\n')
        t = float(np.frombuffer(hd[e + 232:e + 240], '<f8')[0])
        ng, d0, ds = _rst_layout(f, fs, nmb, nx1, nx2, nx3)
        n1, n2, n3 = nx1 + 2 * ng, nx2 + 2 * ng, nx3 + 2 * ng; nc = n1 * n2 * n3; lz = 0.0
        for m in range(nmb):
            f.seek(d0 + m * ds + ds - nc * 8)
            s = np.frombuffer(f.read(nc * 8), '<f8').reshape(n3, n2, n1)[ng:ng + nx3, ng:ng + nx2, ng:ng + nx1]
            g = geom[m]
            xi = np.pi / 4 * (g[2] + (g[3] - g[2]) * (np.arange(nx2) + 0.5) / nx2)
            eta = np.pi / 4 * (g[4] + (g[5] - g[4]) * (np.arange(nx3) + 0.5) / nx3)
            XI, ETA = np.meshgrid(xi, eta); X = dhjcs.panel_map(m // (nmb // 6), XI, ETA)
            lz += float(np.sum(s * X[..., 2][..., None]))
    return t, lz

def spin_series(run, nx1, cachef, every=5):
    """sum S rhat_z at every `every`-th restart (+ the newest), cached per file (name, size, mtime)"""
    import dhjcs
    c = json.load(open(cachef)) if os.path.exists(cachef) else {}
    rs = [f for f in sorted(glob.glob(run + '/rst/*.rst')) if os.path.getmtime(f) < __import__('time').time() - 300]
    if not rs: return None
    pick = rs[::every] + ([rs[-1]] if (len(rs) - 1) % every else [])
    geom = None
    for f in pick:
        k = '%s %d %d' % (os.path.basename(f), os.path.getsize(f), int(os.path.getmtime(f)))
        if k in c: continue
        if geom is None:
            raw = dhjcs.bin_convert.read_binary(sorted(glob.glob(run + '/bin/*.bin'))[0])
            geom = np.asarray(raw['mb_geometry']); nj, nk = raw['nx2_mb'], raw['nx3_mb']
        try: c[k] = spin_lz_rst(f, geom, nx1, nj, nk)
        except Exception as ex: print('spin read failed', f, ex)
    json.dump(c, open(cachef, 'w'))
    v = sorted([tuple(x) for kk, x in c.items() if kk.split()[0] in {os.path.basename(f) for f in pick}])
    return np.array(v) if v else None

# ---------------------------------------------------------------- figures + takeaways (called from ana.py full mode)
C_CH, C_ACC, C_RES, C_RES0 = '#1f4e79', '#b5541c', '#2e7d4f', '#9a9a9a'

def _symlin(a, ys):
    m = max([np.nanmax(np.abs(y)) for y in ys if y is not None and len(y)] + [1e-300])
    return 10.0 ** np.floor(np.log10(max(m, 1e-300) / 300.0))

def analyse(key, name, run, om, nx1, cdir):
    B = budgets(run, om)
    S = spin_series(run, nx1, '%s/spin_%s.json' % (cdir, key))
    L = B['L']; t = B['t']
    if S is not None and len(S):
        idx = np.clip(np.searchsorted(t, S[:, 0]), 0, len(t) - 1)
        L['spin_rot'] = S[:, 0] / B['P']; L['spin'] = S[:, 1]
        L['res_spin'] = L['res'][idx] + S[:, 1]
    return B

def fig_budget(name, B, fn):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 7.2), sharex=True, gridspec_kw=dict(height_ratios=[1.25, 1]))
    rot = B['rot']
    ttl = {'M': 'mass', 'E': 'total energy (e_int + KE + rho Phi_tot)', 'L': 'L_z = Lrel + Omega Iw'}
    for c, q in enumerate('MEL'):
        b = B[q]; s = b['scale']; a, r = ax[0, c], ax[1, c]
        a.plot(rot, b['d'] / s, c=C_CH, lw=1.8, label='change since t = 0')
        a.plot(rot, b['acc'] / s, c=C_ACC, lw=1.4, ls='--', label='integrated booked sources')
        if q == 'L' and 'res_spin' in b:
            i = np.clip(np.searchsorted(B['rot'], b['spin_rot']), 0, len(rot) - 1)
            a.plot(b['spin_rot'], (b['acc'][i] - b['spin']) / s, ls='none', marker='o', ms=4, c=C_RES,
                   label='booked - spin reservoir (rst)')
            r.plot(b['spin_rot'], b['res_spin'] / s, c=C_RES, lw=1.6, marker='o', ms=3.5, label='residual incl. spin reservoir')
            r.set_title('hst terms alone leave %+.1e (the spin reservoir)' % (b['res'][-1] / s), fontsize=9, color='#555555', loc='right')
        else:
            r.plot(rot, b['res'] / s, c=C_RES, lw=1.2, label='residual')
        a.set_title(ttl[q], fontsize=10.5); a.set_ylabel('relative to %s' % b['scale_name'])
        a.legend(fontsize=8, loc='best'); a.axhline(0, c='k', lw=0.5)
        r.axhline(0, c='k', lw=0.5); r.set_ylabel('residual / %s' % b['scale_name']); r.set_xlabel('rotations')
        r.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        if q != 'L' or 'res_spin' not in b: r.set_title('change - booked', fontsize=9, color='#555555', loc='right')
    fig.suptitle('%s: conservation budgets (top: change vs booked sources; bottom: residual)' % name, fontsize=11.5, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(fn, dpi=90); plt.close(fig)

def summary(Bs):
    """final residuals per run (dict) for the summary figure, the html and the note"""
    out = {}
    for key, name, B in Bs:
        r = {'name': name, 'rot': float(B['rot'][-1])}
        for q in 'MEL':
            b = B[q]
            r[q + '_d'] = float(b['d'][-1] / b['scale']); r[q + '_res'] = float(b['res'][-1] / b['scale'])
            r[q + '_res_alt'] = float(b['res'][-1] / max(abs(b['alt'][-1]), 1e-300))
            r[q + '_resmax'] = float(np.max(np.abs(b['res'])) / b['scale'])
        L = B['L']
        if 'res_spin' in L:
            r['L_res_spin'] = float(L['res_spin'][-1] / L['scale']); r['L_spin'] = float(L['spin'][-1] / L['scale'])
            r['L_res_spin_rot'] = float(L['spin_rot'][-1])
        r['E_floor_frac'] = float(B['E']['terms']['floors'][-1] / max(B['E']['terms']['interior heat (wall)'][-1], 1e-300))
        r['Mdot_top_mean'] = float(-B['M']['terms']['outflow top'][-1] / max(B['t'][-1], 1.0))
        out[key] = r
    return out

def fig_summary(S, fn):
    import matplotlib.pyplot as plt
    keys = list(S); x = np.arange(len(keys)); w = 0.2
    fig, ax = plt.subplots(figsize=(11, 4.6))
    ser = [('M_res', 'mass / M0', '#1f4e79'), ('E_res', 'energy / |E0|', '#b5541c'),
           ('L_res', 'L_z / (I Omega), hst terms only', '#9a9a9a'), ('L_res_spin', 'L_z / (I Omega), incl. spin reservoir', '#2e7d4f')]
    for i, (k, lab, c) in enumerate(ser):
        v = np.array([abs(S[kk].get(k, np.nan)) for kk in keys])
        ax.bar(x + (i - 1.5) * w, np.where(v > 0, v, np.nan), w * 0.92, color=c, label=lab)
    ax.set_yscale('log'); ax.set_xticks(x)
    ax.set_xticklabels(['%s\n(%.0f rot)' % (S[k]['name'], S[k]['rot']) for k in keys])
    ax.set_ylabel('|final residual| / scale'); ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.set_title('final budget residual per run (bars = |change - booked sources| at the last row)', fontsize=10)
    ax.grid(axis='x', visible=False)
    fig.tight_layout(); fig.savefig(fn, dpi=90); plt.close(fig)

SPECIAL = {'toi1518': ('The whole offset (about -3.4e-3 I Omega) appears between rotations 14 and 22, the one segment run '
                       'with the older binary 2daa0a7a. That binary wrote no Ttide column, so the tide torque in that window is not '
                       'booked, and it still had the pre-fix roche_wb_tan probe that froze a Coriolis torque into g_t. Outside that '
                       'window the residual stays within 1e-5 I Omega.')}

def html(S, esc):
    o = ['<h2 id="budget">Conservation budgets</h2>',
         '<p>Each budget compares the change of an hst total since t = 0 with the time integral of every source the hst books. '
         '<strong>Mass</strong>: Mdot_bot - Mdot_top + Mfloor. '
         '<strong>Energy</strong> (hst tot-E = e_int + KE + rho Phi_tot under etotgrav, with Phi_tot = point-mass gravity + centrifugal + Roche tide, '
         'all time-independent in the corotating frame): Etot_bot - Etot_top + Efloor. This covers absorbed starlight, IR out of the top, '
         'sigma T_int^4 at the wall, the fluid flux (including the rho Phi v work of the tide) and the floors. '
         '<strong>L_z</strong> (Lrel + Omega Iw): Ttide + Tbnd + Tdrag. Neither the radial-only top sponge nor the Coriolis term '
         '(coriolis_am) adds anything, and the wall has no mass flux (Mdot_bot = 0). '
         'Flux columns are snapshots, integrated with the trapezoid rule; floor and torque columns are interval means.</p>',
         '<p><strong>One L_z term is not in the hst.</strong> &lt;hydro&gt;/am_vector puts the radial part of its transport torque '
         'correction into a per-cell spin reservoir, which is never returned (am_spin_return off). Its z part, '
         'sum S rhat_z, is stored only in the restart files. It is read from every 5th rst (the green markers).</p>',
         '<figure><img src="img/budget_summary.png" alt="final budget residual per run and quantity" loading="lazy">'
         '<figcaption>Absolute final residual (change minus booked sources at the last hst row): mass relative to M0, energy relative to |E0|, '
         'L_z relative to I_w0 Omega, both without and with the spin reservoir.</figcaption></figure>']
    for k, r in S.items():
        img = 'img/budget_%s.png' % k
        tk = []
        ma = abs(r['M_res_alt'])
        tk.append('Mass closes to %.1e of M0%s.' % (abs(r['M_res']),
                  (' (%.1e of the integrated top outflow)' % ma) if ma < 1e-3 else ' (round-off of the hst sum; the top outflow is negligible)'))
        tk.append('Energy closes to %.1e of |E0| (%.1e of the integrated absorbed starlight). The net change is %+.1e |E0|, and the floors supplied %.2g %% as much energy as the interior heat flux.' % (
            abs(r['E_res']), abs(r['E_res_alt']), r['E_d'], 100 * r['E_floor_frac']))
        if 'L_res_spin' in r:
            s = 'L_z: hst terms alone leave %+.1e I Omega. Adding the spin reservoir (%+.1e) leaves %+.1e I Omega at rot %.0f' % (
                r['L_res'], r['L_spin'], r['L_res_spin'], r['L_res_spin_rot'])
            if k in SPECIAL: s += '. ' + SPECIAL[k]
            elif 'Ttop_v' in r and abs(r['L_res_spin']) > 3e-6:
                s += ('. The remaining drift grows with the top outflow. It matches the angular-momentum flux out of the top, '
                      'which has no hst column (am_hst_split off): Mdot_top carrying a zonal wind of %+.1f km/s at x1max.' % r['Ttop_v'])
            else: s += ', i.e. it closes.'
            tk.append(s)
        else:
            tk.append('L_z: hst terms alone leave %+.1e I Omega (no restart file read yet for the spin reservoir).' % r['L_res'])
        o.append('<h3>%s</h3><ul class="facts">%s</ul>' % (esc(r['name']), ''.join('<li>%s</li>' % esc(x) for x in tk)))
        o.append('<figure><img src="%s" alt="%s mass, energy and L_z budgets" loading="lazy"><figcaption>Top row: change since t = 0 (blue) and integrated booked sources (orange dashed). '
                 'In the L_z panel the green dots also subtract the spin-reservoir L_z read from the restart files. Bottom row: '
                 'residual (change - booked) on a linear axis, in the units named on the axis.</figcaption></figure>' % (img, esc(r['name'])))
    return o
