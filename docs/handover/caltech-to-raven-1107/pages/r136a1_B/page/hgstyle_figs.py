"""hgstyle_figs.py (r136a1-1104, 11-05 USER): ADDITIONAL line plots for the R136a1 deep page in the He giant global / RSG E
page style (hgglob_page/src/mk_figs.py: dark theme, same rcParams, colours TC, panel titles, r / R_ph,IC on a log axis,
history panels, flux budget with comoving L_rad + EOS enthalpy + kinetic flux). The page's own figures (r136_page.py:
img/timeseries.png, profiles.png, mdot.png) are left as they are; this writes img/hgstyle_{hist,profiles,fluxbudget}.png.

Reductions follow hgglob_page/src/reduce.py (VER 6) on the spherical-polar wedge (shell means = solid-angle weighted over
the wedge, wedge -> star x 4 pi / Omega): rho-weighted rms v_r and |v_t|, T_gas (run's EOS) and T_rad = (E/a)^1/4,
L_rad comoving = lab M1 flux minus the full (E + P_rr) v_r, L_enth with h = (e + P_gas + E + P_rr)/rho minus the mean
mass flux times <h>, L_kin, mean-flow advection, rms d log rho; P_rr = D_rr E from the m1_vet dump at the same time, else E/3.
Plus the R136a1-only quantities (Gamma, P_rad/P, kappa_R, Mdot_r) and the per-column tau_R = 1 radius.
Each dump is reduced once to <cache>/hg_NNNNN.npz (the red_ caches of r136_page.py are not touched).
usage: python3 hgstyle_figs.py PAGEDIR"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r136_page as rp  # noqa: E402   (Table, EOS, assemble, dumps, x1faces, ic_column, r_photo; RUN / CACHE from env)
bc = rp.bc
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import ScalarFormatter, NullFormatter  # noqa: E402

sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/pagefast_lib')
import hstcache as _hc  # noqa: E402   (the HG/RSG pages' exact incremental hst reader, restart filter included)

# ---- HG / RSG page style (hgglob_page/src/mk_figs.py) ----
plt.rcParams.update({'figure.facecolor': '#0b0d10', 'axes.facecolor': '#12161b', 'savefig.facecolor': '#0b0d10',
                     'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
                     'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'axes.titlesize': 11,
                     'font.size': 10, 'legend.fontsize': 8, 'legend.facecolor': '#12161b',
                     'legend.edgecolor': '#3a424c', 'grid.color': '#262c33'})
TC = ['#9aa4ae', '#4cc3ff', '#7dff8a', '#ffb347', '#ff5c7a']
VER = 3          # 3 (11-06): F0 solver-consistent where G > 1e-3 (r136_f0.py)
A_RAD = 7.5657e-15
DAY, RSUN, LSUN, MSUN = 86400., 6.957e10, 3.828e33, 1.989e33
FAC, G, C = rp.FAC, rp.G, rp.C


def smooth(y, n):
    if len(y) < 2 * n or n < 2:
        return y
    c = np.cumsum(np.insert(y, 0, 0.))
    s = (c[n:] - c[:-n]) / n
    return np.concatenate([np.full(n - 1, np.nan), s])


def vet_map():
    m = {}
    for f in glob.glob(rp.RUN + '/bin/r136.m1_vet.*.bin'):
        try:
            m[round(rp.dump_time(f))] = f
        except Exception:
            pass
    return m


def reduce_hg(n, hf, mf, x1f, kap, eos, vmap):
    os.makedirs(rp.CACHE, exist_ok=True)
    cf = os.path.join(rp.CACHE, 'hg_%05d.npz' % n)
    t_h = rp.dump_time(hf)
    fv = vmap.get(round(t_h))
    st = np.array([os.path.getsize(hf), os.path.getsize(mf), os.path.getsize(fv) if fv else 0], float)
    if os.path.exists(cf):
        c = np.load(cf)
        if int(c['ver']) >= VER and np.array_equal(c['stamp'], st):
            return dict(c)
    dh = bc.read_binary(hf); t = float(dh['time']); w = rp.assemble(dh); del dh
    dm = bc.read_binary(mf); m = rp.assemble(dm); del dm
    nx3, nx2, nx1 = w['dens'].shape
    x2f = np.linspace(np.pi / 3, 2 * np.pi / 3, nx2 + 1); x3f = np.linspace(0, np.pi / 3, nx3 + 1)
    r = 0.5 * (x1f[1:] + x1f[:-1]); dr = np.diff(x1f); dvr = (x1f[1:] ** 3 - x1f[:-1] ** 3) / 3
    dom = np.diff(x3f)[:, None] * (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :]
    om = dom.sum()
    sm = lambda a: np.einsum('kji,kj->i', a.astype(np.float64), dom) / om   # noqa: E731
    rho, vr, v2, v3, ei = w['dens'], w['velx'], w['vely'], w['velz'], w['eint']
    del w
    T, Pg = eos.temp_pgas(rho, ei)
    E, F1 = m['m1_e'], m['m1_f1']
    del m
    kr = kap(rho, T)
    res = dict(ver=VER, stamp=st, t=t, n=n, r=r, x1f=x1f)
    rhom = sm(rho); res['rho'] = rhom
    rvf = rho * vr; mrv = sm(rvf); mvr = mrv / rhom
    res['vrrms'] = np.sqrt(np.maximum(sm(rho * vr ** 2) / rhom - mvr ** 2, 0))
    vt2 = v2 ** 2 + v3 ** 2
    res['vtrms'] = np.sqrt(sm(rho * vt2) / rhom)
    res['Tgas'] = sm(T); res['Trad'] = sm((np.maximum(E, 1e-30) / A_RAD) ** 0.25)
    res['E'] = sm(E); Fr = sm(F1); res['Fr'] = Fr
    Lr = 4 * np.pi * r ** 2 * Fr; res['Lr'] = Lr
    res['kap'] = sm(kr)
    res['prad_frac'] = sm((E / 3) / (E / 3 + Pg))
    # Gamma (as r136_page.py): m(r) = m(x1min) + mass of the dump below r
    dmass = np.einsum('kji,kj->i', rho.astype(np.float64), dom) * dvr
    mcum = rp.GM0 / G + FAC * (np.cumsum(dmass) - 0.5 * dmass)
    g = G * mcum / r ** 2
    res['gam_rad_lab'] = sm(kr * F1) / (C * g)
    res['gam_L'] = res['kap'] * rp.LSTAR / (4 * np.pi * C * G * mcum)
    res['mdot'] = 4 * np.pi * r ** 2 * mrv
    res['M_env'] = FAC * dmass.sum()
    # HG reduce.py VER 6 flux split
    hid = (5. / 3.) * ei / rho
    res['Lconv_id'] = 4 * np.pi * r ** 2 * (sm(rvf * hid) - mrv * sm(hid))
    del hid
    hent = ei / rho + Pg / rho
    Lconv_gas = 4 * np.pi * r ** 2 * (sm(rvf * hent) - mrv * sm(hent))
    Lbulk = 4 * np.pi * r ** 2 * mrv * sm(hent)
    vet = False
    if fv:
        try:
            hv = bc.read_binary(fv)
            if abs(float(hv['time']) - t) < 1e-6 * max(abs(t), 1.0):
                hent = (1.0 + rp.assemble(hv)['m1_d11']) * E / rho; vet = True
            del hv
        except Exception as e:
            print('m1_vet unreadable', n, e, flush=True)
    if not vet:
        hent = (4. / 3.) * E / rho
    F0, _ = rp.f0_solver(F1 - rvf * hent, E, kr, rho, x1f, r)  # comoving F0 (solver-consistent, r136_f0.py)
    res['gam_rad'] = sm(kr * F0) / (C * g)
    Lrcom = 4 * np.pi * r ** 2 * sm(F0)
    Lradbulk = 4 * np.pi * r ** 2 * mrv * sm(hent)
    Lradadv = Lr - Lrcom - Lradbulk                            # keeps L_rad + L_enth + L_kin = the lab total
    del hent, F0
    res['Lconv_gas'] = Lconv_gas; res['Lradadv'] = Lradadv
    res['Lconv'] = Lconv_gas + Lradadv
    res['Lrcom'] = Lrcom
    ke = 0.5 * (vr ** 2 + vt2)
    res['Lkin'] = 4 * np.pi * r ** 2 * (sm(rvf * ke) - mrv * sm(ke))
    res['Lbulk'] = Lbulk + Lradbulk + 4 * np.pi * r ** 2 * mrv * sm(ke)
    res['vet'] = vet
    del ke, rvf, vt2
    lrho = np.log10(rho)
    res['dlrho_rms'] = np.sqrt(sm((lrho - sm(lrho)[None, None, :].astype(np.float32)) ** 2))
    del lrho
    # tau_R from the top inward: shell mean (page R_ph) and per column (spread)
    dtau = (kr * rho).astype(np.float64) * dr[None, None, :]
    tau = np.flip(np.cumsum(np.flip(dtau, axis=2), axis=2), axis=2)
    del dtau
    res['tau'] = sm(tau)
    nab = (tau >= 1.0).sum(axis=2)                      # tau decreases outward: last index with tau >= 1 is nab - 1
    i = np.clip(nab - 1, 0, nx1 - 2)
    ta = np.take_along_axis(tau, i[..., None], 2)[..., 0]; tb = np.take_along_axis(tau, (i + 1)[..., None], 2)[..., 0]
    la, lb = np.log(np.maximum(ta, 1e-300)), np.log(np.maximum(tb, 1e-300))
    rc = r[i] + (r[i + 1] - r[i]) * np.clip((0 - la) / np.where(lb != la, lb - la, -1), 0, 1)
    rc = np.where(nab == 0, np.nan, np.where(nab >= nx1, r[-1], rc))
    res['Rph_p2'], res['Rph_med'], res['Rph_p98'] = (float(np.nanpercentile(rc, q)) for q in (2, 50, 98))
    del tau
    np.savez(cf, **res)
    print('hgstyle reduced dump %05d t = %.3f d (vet %s)' % (n, t / DAY, vet), flush=True)
    return res


def logx(a):
    a.set_xscale('log')
    for ax_ in (a.xaxis,):
        ax_.set_major_formatter(ScalarFormatter()); ax_.set_minor_formatter(NullFormatter())
    a.set_xticks([0.8, 1, 1.5, 2, 3, 4])


def main():
    page = sys.argv[1]
    out = lambda s: os.path.join(page, 'img', 'hgstyle_%s.png' % s)  # noqa: E731
    x1f = rp.x1faces(); kap = rp.Table(rp.OPAC); eos = rp.EOS(rp.EOSF); ic = rp.ic_column()
    ic['Pg'] = eos.temp_pgas(ic['rho'], ic['eint'])[1].astype(np.float64)
    vmap = vet_map()
    R = [reduce_hg(n, hf, mf, x1f, kap, eos, vmap) for n, hf, mf in rp.dumps()]
    R.sort(key=lambda x: float(x['t']))
    # R_ph of the IC column: tau_R = 1 integrated inward from the column top (same opacity table)
    icr = ic['r']; dri = np.gradient(icr)
    taui = np.flip(np.cumsum(np.flip(kap(ic['rho'], ic['T']).astype(np.float64) * ic['rho'] * dri)))
    RPH0 = rp.r_photo(icr, taui)
    RPHS = '%.1f R_sun' % (RPH0 / RSUN)
    XIN, XOUT = x1f[0] / RPH0, x1f[-1] / RPH0
    H = _hc.readhst(rp.RUN + '/r136.user.hst')
    td = H['time'] / DAY
    Lin = FAC * H['L_in']
    tr = np.array([float(d['t']) for d in R]) / DAY

    # ---------------- history (fig_hist) ----------------
    fig, ax = plt.subplots(2, 3, figsize=(15, 7.8), constrained_layout=True)
    n = max(1, min(100, len(td) // 60))
    a = ax[0, 0]
    for k, c, lab in (('L_top', '#ffb347', 'L_top'), ('L_mid', '#c77dff', 'L_mid'), ('L_bot', '#4cc3ff', 'L_bot'),
                      ('L_MLT', '#7dff8a', 'L_MLT (frozen scaffold)')):
        a.plot(td, FAC * H[k] / Lin, color=c, lw=1.0, label=lab + ' / L_in')
    a.axhline(1, color='#a3adb8', lw=0.6, ls=':')
    a.set_ylabel('luminosity / L_in'); a.legend(loc='best')
    a.set_title('luminosities (L_in = %.3g L_sun imposed at r_in)' % (Lin[-1] / LSUN))
    a = ax[0, 1]
    a.semilogy(td, FAC * np.maximum(H['KE_int'], 1e30), color='#ff5c7a', lw=1.2, label='KE_int')
    a.semilogy(td, FAC * np.maximum(H['KEr_int'], 1e30), color='#c77dff', lw=1.0, label='KE_r,int (radial)')
    a.set_ylabel('erg'); a.legend(); a.set_title('kinetic energy inside the star')
    a = ax[0, 2]
    E0 = H['Etot'][0]
    a.plot(td, (H['Etot'] - E0) / abs(E0), color='#e8ecf1', lw=1.2, label='(E_tot − E_tot,0)/|E_tot,0| (user hst)')
    ct = np.concatenate([[0], np.cumsum(0.5 * (H['L_in'][1:] + H['L_in'][:-1] - H['L_top'][1:] - H['L_top'][:-1]) * np.diff(H['time']))])
    a.plot(td, ct / abs(E0), color='#ffb347', lw=1.0, ls='--', label='∫(L_in − L_top) dt / |E_tot,0|')
    M0 = H['M_tot'][0]
    a.plot(td, (H['M_tot'] - M0) / M0, color='#4cc3ff', lw=1.0, label='(M_tot − M_0)/M_0')
    a.set_title('energy and mass budget'); a.legend(loc='best'); a.set_ylabel('relative change')
    a = ax[1, 0]
    a.plot(td, H['dt'], color='#7dff8a', lw=1.0)
    a.set_yscale('log'); a.set_ylabel('dt [s]'); a.set_title('time step (log)')
    a2 = a.twinx()
    a2.plot(td, smooth(H['Picard'], n), color='#e8ecf1', lw=0.8, alpha=0.8)
    a2.set_ylabel('Picard per cycle (running mean, white)')
    a = ax[1, 1]
    rph = np.array([rp.r_photo(d['r'], d['tau']) for d in R])
    if len(tr):
        a.fill_between(tr, [float(d['Rph_p2']) / RPH0 for d in R], [float(d['Rph_p98']) / RPH0 for d in R],
                       color='#4cc3ff', alpha=0.25, lw=0, label='per-column τ_R = 1: 2–98 %')
        a.plot(tr, [float(d['Rph_med']) / RPH0 for d in R], color='#4cc3ff', lw=0.8, ls='--', label='per-column median')
        a.plot(tr, rph / RPH0, 'o-', ms=2.5, color='#4cc3ff', lw=1.2, label='R_ph (shell-mean τ_R = 1)')
    a.axhline(1, color='#a3adb8', lw=0.6, ls=':')
    a.set_ylabel('R_ph / R_ph,IC (%s)' % RPHS); a.legend(loc='upper left'); a.set_title('photosphere (from the 0.25 d dumps)')
    a2 = a.twinx()
    if len(tr):
        Lt = np.interp(tr * DAY, H['time'], FAC * H['L_top'])
        a2.plot(tr, (Lt / (4 * np.pi * rp.SIG * rph ** 2)) ** 0.25, color='#ffb347', lw=0.9, ls='--')
    a2.set_ylabel('T_eff [K] from L_top, R_ph (dashed)', color='#ffb347')
    a = ax[1, 2]
    a.semilogy(td, FAC * np.maximum(np.abs(H['Mdot_top']), 1e8), color='#66d9c6', lw=1.0, label='|Mdot_top| [g/s]')
    a.semilogy(td, FAC * np.maximum(np.abs(H['P_sponge']), 1e20), color='#c77dff', lw=1.0, label='|P_sponge| [erg/s]')
    a.semilogy(td, FAC * np.maximum(np.abs(H['P_esrc']), 1e20), color='#ffb347', lw=1.0, label='|P_esrc| [erg/s]')
    a.set_title('outer boundary / sources'); a.legend()
    for a in ax.flat:
        a.set_xlabel('t [d]'); a.grid(True, lw=0.5)
    fig.savefig(out('hist'), dpi=90)
    plt.close(fig)

    # ---------------- radial profiles (fig_profiles + R136a1 row) ----------------
    icx = ic['r'] / RPH0
    pick = rp.pick_times([None] * len(R))
    wm = np.interp([float(R[p]['t']) for p in pick], H['time'], H['w_mlt'])
    fig, ax = plt.subplots(3, 3, figsize=(15, 12), constrained_layout=True)
    LN = float(Lin[-1])
    for c, (p, wmi) in enumerate(zip(pick, wm)):
        d = R[p]; x = d['r'] / RPH0; col = TC[c % len(TC)]
        lab = 't = %.2f d' % tr[p]
        ax[0, 0].loglog(x, d['rho'], color=col, lw=1, label=lab)
        ax[0, 1].loglog(x, d['Tgas'], color=col, lw=1.1, label=lab)
        ax[0, 1].loglog(x, d['Trad'], color=col, lw=0.6, ls='--')
        ax[0, 2].semilogx(x, d['vrrms'] / 1e5, color=col, lw=1.1, label=lab + ' v_r')
        ax[0, 2].semilogx(x, d['vtrms'] / 1e5, color=col, lw=0.8, ls='--')
        ax[1, 0].semilogx(x, d['Lr'] / LN, color=col, lw=0.5, alpha=0.3)
        ax[1, 0].semilogx(x, d['Lrcom'] / LN, color=col, lw=1.1, label=lab)
        ax[1, 1].semilogx(x, d['Lconv'] / LN, color=col, lw=1.1, label=lab)
        ax[1, 1].semilogx(x, d['Lkin'] / LN, color=col, lw=0.8, ls='--')
        if wmi > 1e-3:
            ax[1, 1].semilogx(icx, wmi * ic['fmlt'], color=col, lw=0.7, ls=':')
        ax[1, 2].semilogx(x, d['dlrho_rms'], color=col, lw=1.1, label=lab)
        ax[2, 0].semilogx(x, d['gam_rad_lab'], color=col, lw=0.5, alpha=0.3)
        ax[2, 0].semilogx(x, d['gam_rad'], color=col, lw=1.1, label=lab)
        ax[2, 0].semilogx(x, d['gam_L'], color=col, lw=0.7, ls=':')
        ax[2, 1].semilogx(x, d['prad_frac'], color=col, lw=1.1, label=lab)
        ax[2, 2].semilogx(x, d['mdot'] / d['M_env'] * DAY, color=col, lw=1.1, label=lab)
    ax[0, 0].loglog(icx, ic['rho'], color='w', lw=0.9, ls=':', label='IC column')
    ax[0, 1].loglog(icx, ic['T'], color='w', lw=0.9, ls=':', label='IC T_gas')
    ax[1, 1].semilogx(icx, ic['fmlt'], color='w', lw=0.9, ls=':', label='IC scaffold F_MLT/F (× w_mlt(t): coloured dotted)')
    ax[2, 1].semilogx(icx, ic['E'] / 3 / (ic['E'] / 3 + ic['Pg']), color='w', lw=0.9, ls=':', label='IC column')
    ax[0, 0].set_ylabel('⟨ρ⟩ [g cm⁻³]'); ax[0, 0].set_title('density (shell mean)')
    ax[0, 1].set_ylabel('⟨T⟩ [K]'); ax[0, 1].set_title('gas temperature T_gas (solid), T_rad (dashed); shell means')
    ax[0, 2].set_ylabel('km/s'); ax[0, 2].set_title('ρ-weighted rms v_r (solid), rms |v_t| (dashed)')
    ax[1, 0].set_ylabel('L_rad(r) / L_in'); ax[1, 0].set_title('comoving L_rad 4πr²⟨F_r − (E+P_rr)v_r⟩ (faint: lab frame)', fontsize=10); ax[1, 0].set_ylim(-0.2, 1.6)
    ax[1, 1].set_ylabel('L_conv,res(r) / L_in'); ax[1, 1].set_title('resolved enthalpy flux, gas+rad h (solid); kinetic-energy flux (dashed)', fontsize=10); ax[1, 1].set_ylim(-0.3, 1.2)
    ax[1, 2].set_yscale('log'); ax[1, 2].set_ylabel('rms δ log ρ'); ax[1, 2].set_title('lateral density fluctuation per shell')
    gmax = max(2.0, min(6.0, 1.1 * max(float(np.nanmax(R[p]['gam_rad'])) for p in pick))) if pick else 2.0
    ax[2, 0].axhline(1, color='#a3adb8', lw=0.6, ls='--'); ax[2, 0].set_ylim(0, gmax)
    ax[2, 0].set_ylabel('Γ'); ax[2, 0].set_title('Γ = ⟨κ_R F_0⟩/(c g), comoving (solid; faint: lab), κ_R L/(4πcGm) (dotted)', fontsize=10)
    ax[2, 1].set_ylim(0, 1); ax[2, 1].set_ylabel('P_rad / (P_rad + P_gas)'); ax[2, 1].set_title('radiation pressure fraction (P_rad = E/3)')
    ax[2, 2].axhline(0, color='#a3adb8', lw=0.5); ax[2, 2].set_yscale('symlog', linthresh=1e-4)
    ax[2, 2].set_ylabel('Ṁ_r / M_env [d⁻¹]  (> 0 outward)')
    ax[2, 2].set_title('mass flux 4πr²⟨ρv_r⟩ / M_env (M_env = %.3g M_sun)' % (float(R[0]['M_env']) / MSUN) if R else 'mass flux')
    for a in ax.flat:
        a.axvline(1, color='#a3adb8', lw=0.5, ls=':')
        a.axvline(XIN, color='#a3adb8', lw=0.5, ls='--')
        a.set_xlim(XIN * 0.98, XOUT); logx(a)
        a.set_xlabel('r / R_ph,IC (R_ph,IC = %s; dashed: x1min)' % RPHS); a.grid(True, lw=0.5); a.legend(fontsize=7)
    fig.savefig(out('profiles'), dpi=100)
    plt.close(fig)

    # ---------------- energy-flux budget (fig_fluxbudget) ----------------
    have = [p for p in range(len(R)) if bool(R[p]['vet'])]
    if len(have) >= 2:
        nw = max(1, min(10, len(have) // 4))
        wins = [have[len(have) // 2 - nw // 2: len(have) // 2 - nw // 2 + nw], have[-nw:]]
        fig, ax = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True, sharey=True)
        for a, wl in zip(ax, wins):
            x = R[wl[-1]]['r'] / RPH0
            av = lambda k: np.mean([np.asarray(R[p][k]) for p in wl], axis=0) / LN  # noqa: E731
            lr, le, lk, lb = av('Lrcom'), av('Lconv'), av('Lkin'), av('Lbulk')
            a.semilogx(x, lr, color='#4cc3ff', lw=1.1, label='L_rad, comoving (lab M1 flux − (E + P_rr) v_r)')
            a.semilogx(x, le, color='#ffb347', lw=1.3, label='L_enth, h = (e + P_gas + E + P_rr)/ρ')
            a.semilogx(x, lk, color='#ff5c7a', lw=1.3, label='L_kin')
            a.semilogx(x, av('Lconv_gas'), color='#ffb347', lw=0.8, ls='--', label='L_enth, gas part only (EOS)')
            a.semilogx(x, av('Lconv_id'), color='#ffb347', lw=0.8, ls=':', label='L_enth, ideal-gas h')
            a.semilogx(x, lb, color='#9aa4ae', lw=1.0, ls='--', label='mean-flow advection 4πr²⟨ρv_r⟩⟨h + v²/2⟩')
            tot = lr + le + lk; tl = 'L_rad + L_enth + L_kin'
            wmean = float(np.mean(np.interp([float(R[p]['t']) for p in wl], H['time'], H['w_mlt'])))
            if wmean > 1e-3:
                lm = wmean * np.interp(R[wl[-1]]['r'], ic['r'], ic['fmlt'])
                a.semilogx(x, lm, color='#7dff8a', lw=1.0, label='L_MLT scaffold (w_mlt = %.2f)' % wmean)
                tot = tot + lm; tl += ' + L_MLT'
            a.semilogx(x, tot, color='#e8ecf1', lw=1.6, label=tl)
            a.semilogx(x, tot + lb, color='#e8ecf1', lw=0.9, ls=':', label='… + mean-flow advection (= lab-frame total)')
            a.axhline(1, color='#a3adb8', lw=0.8, ls='--'); a.axhline(0, color='#a3adb8', lw=0.5)
            a.axvline(1, color='#a3adb8', lw=0.5, ls=':')
            a.set_title(('mean of %d dumps, t = %.2f–%.2f d' % (len(wl), tr[wl[0]], tr[wl[-1]])) if len(wl) > 1
                        else 'dump at t = %.2f d' % tr[wl[0]])
            a.set_xlim(XIN * 0.98, XOUT); logx(a)
            a.set_xlabel('r / R_ph,IC'); a.grid(True, lw=0.5); a.legend(fontsize=8)
        ax[0].set_ylabel('luminosity / L_in'); ax[0].set_ylim(-0.6, 1.8)
        fig.savefig(out('fluxbudget'), dpi=100)
        plt.close(fig)
    print('hgstyle figs done: %d dumps (%d with m1_vet), R_ph,IC = %s' % (len(R), len(have), RPHS))


if __name__ == '__main__':
    main()
