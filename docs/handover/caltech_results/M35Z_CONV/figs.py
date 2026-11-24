#!/usr/bin/env python3
"""figures + summary.txt for M35Z_CONV.  usage: figs.py REDDIR OUTDIR"""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import red as R
RD, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
G = R.GM; RS = 6.957e10; LS = 6.618281e38
Z = [dict(np.load(f)) for f in sorted(glob.glob(RD + '/m35_*.npz'))]
Z.sort(key=lambda z: float(z['t']))
t = np.array([float(z['t']) for z in Z]); r = Z[-1]['r']; x = r / RS
iE = int(np.argmin(abs(t - 19800.))); iL = len(Z) - 1
TXT = []
def p(*a):
    s = ' '.join(str(v) for v in a); print(s, flush=True); TXT.append(s)

def state(z):
    """shell-mean structure: nabla, nabla_ad, N2, from mean T, rho, Pg, E (P = Pg + E/3; ideal gas gamma 5/3 + radiation)"""
    P = z['Pg'] + z['Er'] / 3.; b = z['Pg'] / P
    G1 = b + (4 - 3 * b) ** 2 * (2 / 3.) / (b + 8 * (1 - b))
    nad = (8 - 6 * b) / (32 - 24 * b - 3 * b * b)
    lP, lT, lr_ = np.log(P), np.log(z['T']), np.log(z['rho'])
    nab = np.gradient(lT, r) / np.gradient(lP, r)
    Trad = (z['Er'] / R.A) ** 0.25; nabr = np.gradient(np.log(Trad), r) / np.gradient(lP, r)
    g = G / r ** 2
    N2 = -g * (np.gradient(lr_, r) - np.gradient(lP, r) / G1)
    Hp = -1 / np.gradient(lP, r)
    cs = np.sqrt(G1 * P / z['rho'])
    Gam = z['kapR'] * (z['Lrad_com'] / (4 * np.pi * r ** 2)) / (R.C * g)
    chi = z['kapR'] * z['rho']
    # radiative diffusivity of the gas+radiation mixture [cm^2/s], estimate: K = 4 a c T^3/(3 chi); cp = d(h)/dT at const P
    cv = 1.5 * z['Pg'] / z['rho'] / z['T'] + 4 * z['Er'] / z['rho'] / z['T']
    dlt = (4 - 3 * b) / b; cp = cv + (P / (z['rho'] * z['T'])) * (4 - 3 * b) ** 2 / b
    Krad = 4 * R.A * R.C * z['T'] ** 3 / (3 * chi); kt = Krad / (z['rho'] * cp)
    return dict(P=P, b=b, G1=G1, nad=nad, nab=nab, nabr=nabr, N2=N2, Hp=Hp, cs=cs, Gam=Gam, kt=kt, g=g, dlt=dlt, Trad=Trad)
SE, SL = state(Z[iE]), state(Z[iL])
zl = Z[iL]
# iron peak: shell in T = 1.2e5..2.6e5 K with the largest kappa_R (mean), on the last dump
sel = (zl['T'] > 1.2e5) & (zl['T'] < 2.6e5)
iFe = int(np.where(sel)[0][np.argmax(zl['kapR'][sel])])
tau = zl['tau_sh']
def at_tau(lev): return int(len(tau) - 1 - np.argmax(tau[::-1] >= lev))
i1, i23 = at_tau(1.0), at_tau(2 / 3.)
def rtau(z, lev):                                  # shell tau at inner faces (log-linear in r)
    tt = z['tau_sh']; k = len(tt) - 1 - np.argmax(tt[::-1] >= lev); rf = z['rf']
    return (rf[k] + (rf[k + 1] - rf[k]) * (np.log(tt[k]) - np.log(lev)) / (np.log(tt[k]) - np.log(tt[k + 1]))) / RS

p('# M35ZAMS convection diagnostic (T2 run, prod_rc4v_1124/m35zams, binary athena_ia123rc4v_262affd6f). NOTE 2026-11-24zzzo.')
p('# dumps %d, t = %.0f .. %.0f s (%.3f d); early dump t = %.0f s, last dump t = %.0f s; cgs unless stated' % (len(Z), t[0], t[-1], t[-1] / 86400, t[iE], t[iL]))
p('# all quantities are shell (solid-angle) statistics of the wedge; L scaled to 4 pi; rms are rho-weighted. ESTIMATES: nabla, N2 (from shell-mean T, rho, Pg, E; P = Pg + E/3), comoving L_rad and convective fluxes (cell F1, E; P_rr = D11 E on the vet_vet dumps t = 0, 43200, 77110 s, E/3 elsewhere), diffusivity, growth/damping times')
p('# check E/(a T_gas^4) deep (i = 50) at last dump: %.4f (1 = LTE)' % zl['E_med_over_aT4'])
p('r_top = %.4f R_sun (tau integrated inward from the top face). shell-mean tau = 1 at r = %.4f (early) / %.4f (last) R_sun, tau = 2/3 at %.4f / %.4f; per-column tau = 1 radius, last dump, mean/min/max: %s R_sun' % (zl['rf'][-1] / RS, rtau(Z[iE], 1.), rtau(zl, 1.), rtau(Z[iE], 2 / 3.), rtau(zl, 2 / 3.), np.round(zl['rc_tau1'] / RS, 4)))
p('iron peak: i = %d r = %.4f R_sun, T = %.4e K, kappa_R = %.4f cm2/g (shell T range 1.2-2.6e5 K, kappa maximum)' % (iFe, x[iFe], zl['T'][iFe], zl['kapR'][iFe]))

def row(nm, i, z, S):
    v = lambda k: float(z[k][i])
    return '%-9s %-7.4f %.3e %.3e %.3e %9.4f %8.3f %7.2f %8.4f %8.4f %7.2f %8.3f %9.3f %9.3f %9.4f %9.4f %.2e' % (
        nm, x[i], v('rho'), v('T'), float(z['tau_sh'][i]), S['nab'][i] - S['nad'][i], S['N2'][i] / (S['g'][i] / S['Hp'][i]), 1e-5 * v('vr_rms_fl'), 1e-5 * v('vt_rms'), 1e-5 * v('vp_rms'), v('mach_t_rms') * 1e3, v('mach_g_max'), v('kapR'), S['Gam'][i], v('Lrad_com') / LS, v('Lenth') / LS, v('Lkin') / LS)
hdr = '%-9s %-7s %-9s %-9s %-9s %9s %8s %7s %8s %8s %7s %8s %9s %9s %9s %9s %s' % ('mark', 'r/Rsun', 'rho', 'T[K]', 'tau', 'nab-nad', 'N2/(g/H)', 'vr_fl', 'vth', 'vph', '1e3Mtot', 'Mg_max', 'kappa_R', 'Gamma', 'Lradc/L', 'Lenth/L', 'Lkin/L')
p('\nTable: rms velocities in km/s (vr_fl = about the shell mass-weighted mean), N2/(g/H_P) dimensionless (<0 = unstable to adiabatic convection by the shell-mean gradients)')
for lab, z, S in (('early t=%.0f' % t[iE], Z[iE], SE), ('last t=%.0f' % t[iL], zl, SL)):
    p('\n' + lab); p(hdr)
    for nm, i in (('Fe peak', iFe), ('r(tau1)', i1), ('r=7.83', int(np.argmin(abs(x - 7.83)))), ('r=7.6', int(np.argmin(abs(x - 7.6)))), ('r=7.4', int(np.argmin(abs(x - 7.4)))), ('r=7.0', int(np.argmin(abs(x - 7.0))))):
        p(row(nm, i, z, S))
# tau<1 excluded
p('\nSame rms with cells of column tau < 1 excluded (rho-weighted, last dump; NaN = no cells)')
p('%-9s %-7s %8s %8s %8s %8s %8s' % ('mark', 'r/Rsun', 'vr_fl', 'vth', 'vph', '1e3Mtot', 'dT_rms'))
for nm, i in (('Fe peak', iFe), ('r=7.83', int(np.argmin(abs(x - 7.83)))), ('r(tau1)', i1), ('r(tau1)+1', i1 + 1), ('r(tau1)+2', i1 + 2), ('r=7.4', int(np.argmin(abs(x - 7.4))))):
    z = zl; p('%-9s %-7.4f %8.4f %8.4f %8.4f %8.3f %8.2e' % (nm, x[i], 1e-5 * z['vr_rms_fl_m1'][i], 1e-5 * z['vt_rms_m1'][i], 1e-5 * z['vp_rms_m1'][i], 1e3 * z['mach_t_rms_m1'][i], z['dT_rms_m1'][i]))

for lab, z, S in (('early', Z[iE], SE), ('last', zl, SL)):
    u = (S['N2'] < 0) & (z['tau_sh'] > 1); k = np.where(u)[0]
    p('%s: shells with N^2 < 0 at tau > 1: %d cells, r = %s R_sun (T %s K)' % (lab, len(k), '%.4f-%.4f' % (x[k[0]], x[k[-1]]) if len(k) else '-', '%.2e-%.2e' % (z['T'][k[0]], z['T'][k[-1]]) if len(k) else '-'))
    w = z['tau_sh'] > 1; q = z['Lenth'] + z['Lkin']; a = int(np.argmax(np.where(w, q, -1e99)))
    p('%s: max of (L_enth + L_KE)/L at tau > 1 = %.3e at r = %.4f; max |v_r mean-subtracted rms| = %.2f km/s; min L_rad,com/L at tau > 1 = %.5f (cell estimate)' % (lab, q[a] / LS, x[a], 1e-5 * np.nanmax(z['vr_rms_fl'][w]), np.min(z['Lrad_com'][w]) / LS))
p('thin gas (tau << 1): cell-centred L_rad estimate = %.4f L at i = 350, while user.hst L_top/L_in = 1.0002; the two disagree by 4.5 %% and the cause is not checked, so L_rad at tau < 1 is not usable here' % (zl['Lrad_com'][350] / LS))
# ------------ KE by radial shell against time (4 pi scaled)
edges = [0, 40, 90, 150, 200, 250, 290, 320, 345, 360, 372, 384]
bands = [(a, b) for a, b in zip(edges[:-1], edges[1:])]
KF = np.array([[z['KEr_fl'][a:b].sum() + z['KElat_fl'][a:b].sum() for a, b in bands] for z in Z])
KM = np.array([[z['KEr_mean'][a:b].sum() for a, b in bands] for z in Z])
KL = np.array([[z['KElat'][a:b].sum() for a, b in bands] for z in Z])
KTOT = np.array([(z['KEr_fl'] + z['KElat_fl']).sum() for z in Z]); KMT = np.array([z['KEr_mean'].sum() for z in Z]); KLM = np.array([(z['KElat'] - z['KElat_fl']).sum() for z in Z]); KRF = np.array([z['KEr_fl'].sum() for z in Z])
p('\nKE by radial shell (4 pi scaled, erg), t in s; KE_fl = radial + lateral fluctuations about the shell means (shell-mean lateral flow NOT included, listed as KE_latmean); KE_mean = mean radial flow of the shell')
p('e-fold fit: ln KE_fl linear in t over the last %d dumps from t > 2e4 s, s (n. = not monotone), 1/rate in s' % len(t[t > 2e4]))
p('%-18s %12s %12s %12s %10s %10s %8s %s' % ('shell r/Rsun', 'KE_fl early', 'KE_fl last', 'KE_mean last', 'ratio', 'e-fold[s]', 'R^2', 'trend'))
m = t > 2e4
for j, (a, b) in enumerate(bands):
    y = np.log(np.maximum(KF[:, j], 1e-60)); A_ = np.polyfit(t[m], y[m], 1); pred = np.polyval(A_, t[m]); ss = 1 - ((y[m] - pred) ** 2).sum() / max(((y[m] - y[m].mean()) ** 2).sum(), 1e-99)
    tail = y[m][-5:].mean() - y[m][-10:-5].mean()
    p('%-18s %12.3e %12.3e %12.3e %10.2f %10.2e %8.3f %s' % ('%.3f-%.3f' % (x[a], x[min(b, 383)]), KF[iE, j], KF[iL, j], KM[iL, j], KF[iL, j] / max(KF[iE, j], 1e-99), 1 / A_[0] if A_[0] != 0 else np.inf, ss, ('last-10 rise' if tail > 0.02 else 'last-10 fall' if tail < -0.02 else 'last-10 flat')))
p('total KE_fl early %.3e, last %.3e; mean-flow KE early %.3e last %.3e' % (KTOT[iE], KTOT[iL], KMT[iE], KMT[iL]))
p('KE_latmean (shell-mean lateral flow) early %.3e last %.3e; radial part of KE_fl early %.3e last %.3e' % (KLM[iE], KLM[iL], KRF[iE], KRF[iL]))
p('fl. KE over first 5 dumps: %s' % np.round(KTOT[:5], 2))
h = np.loadtxt('/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/m35zams/m35zams.hydro.hst'); hu = np.loadtxt('/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/m35zams/m35zams.user.hst')
idx = np.where(np.diff(h[:, 0]) <= 0)[0]; h = h[idx[-1] + 1:] if len(idx) else h
p('hydro.hst KE 1/2/3 (wedge, code units erg) at t = %.0f: %.3e %.3e %.3e; at last t = %.0f: %.3e %.3e %.3e' % (h[0, 0], h[0, 7], h[0, 8], h[0, 9], h[-1, 0], h[-1, 7], h[-1, 8], h[-1, 9]))
for lab, lo in (('t>2e4', 2e4), ('t>5e4', 5e4)):
    mm = h[:, 0] > lo; s = [np.polyfit(h[mm, 0], np.log(np.maximum(h[mm, k], 1e-60)), 1)[0] for k in (7, 8, 9)]
    p('hydro.hst KE component log-slopes %s: e-fold times (s) 1-KE %.2e, 2-KE %.2e, 3-KE %.2e (negative = decaying)' % ((lab,) + tuple(1 / v if v != 0 else np.inf for v in s)))
KW = 0.5  # (KE e-fold; velocity e-fold is twice this)
p('\nLocal e-fold times of KE [s] by window (slope of ln KE; velocity e-fold = 2x; growth if > 0)')
p('%-12s %12s %12s %12s %12s' % ('window [ks]', 'Fe band', 'hst 1-KE', 'hst 2-KE', 'hst 3-KE'))
jF = [j for j, (a, b) in enumerate(bands) if a <= iFe < b][0]
for lo, hi in ((20, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 78)):
    md = (t >= lo * 1e3) & (t < hi * 1e3); mh = (h[:, 0] >= lo * 1e3) & (h[:, 0] < hi * 1e3)
    f = lambda a, b: 1 / np.polyfit(a, np.log(b), 1)[0]
    p('%-12s %12.3e %12.3e %12.3e %12.3e' % ('%d-%d' % (lo, hi), f(t[md], KF[md, jF]), f(h[mh, 0], h[mh, 7]), f(h[mh, 0], h[mh, 8]), f(h[mh, 0], h[mh, 9])))
p('wedge check: hydro.hst KE sum x 4pi/Omega = %.4e erg vs dumps KE_fl + KE_latmean + KE_mean = %.4e erg' % ((h[-1, 7] + h[-1, 8] + h[-1, 9]) * 4 * np.pi / (2 * np.sin(0.5 * (1.6028405718615124 - 1.5387520817282807)) * 0.064), KTOT[iL] + KLM[iL] + KMT[iL]))
# ------------ growth/damping at the iron peak
for nm, S, z in (('early', SE, Z[iE]), ('last', SL, zl)):
    i = iFe; n2 = S['N2'][i]
    p('%s iron peak (i=%d): N^2 = %.3e s^-2 -> %s; H_P = %.3e cm; thermal diffusivity %.2e cm2/s, H_P^2/kappa_T = %.2e s; c_s(tot) = %.3e cm/s; c/c_s = %.0f; tau_shell(peak) = %.1f (top->inward)' % (
        nm, i, n2, ('unstable, adiabatic growth time 1/sqrt(-N2) = %.2e s' % (1 / np.sqrt(-n2)) if n2 < 0 else 'stable (shell-mean)'), S['Hp'][i], S['kt'][i], S['Hp'][i] ** 2 / S['kt'][i], S['cs'][i], R.C / S['cs'][i], z['tau_sh'][i]))
open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(TXT) + '\n')
np.savetxt(os.path.join(OUT, 'ke_bands.txt'), np.column_stack([t] + [KF[:, j] for j in range(len(bands))]), header='t[s] KE_fl per band (erg, 4pi): ' + ' '.join('%.3f-%.3f' % (x[a], x[min(b, 383)]) for a, b in bands))
# radial profile table (last dump) every 8th cell
np.savetxt(os.path.join(OUT, 'profile_last.txt'), np.column_stack([x, zl['rho'], zl['T'], zl['tau_sh'], zl['kapR'], SL['nab'] - SL['nad'], SL['N2'], 1e-5 * zl['vr_rms_fl'], 1e-5 * zl['vt_rms'], 1e-5 * zl['vp_rms'], zl['mach_t_rms'], zl['mach_t_max'], zl['Lrad_com'] / LS, zl['Lenth'] / LS, zl['Lkin'] / LS])[::4], header='r/Rsun rho T tau kappa_R nabla-nabla_ad N2 vr_fl[km/s] vth vph Mtot_rms Mtot_max Lrad_com/L Lenth/L Lkin/L (last dump, every 4th cell)', fmt='%.5e')

# ------------ figures
def marks(a, xs=True):
    a.axvline(x[i1], color='#1f77b4', lw=0.8, ls=':'); a.axvline(x[iFe], color='#d62728', lw=0.8, ls='--')
fig, ax = plt.subplots(2, 4, figsize=(19, 8.4), constrained_layout=True)
A = ax.ravel()
for z, S, ls, lab in ((Z[iE], SE, '--', 't=%.0f s' % t[iE]), (zl, SL, '-', 't=%.0f s' % t[iL])):
    A[0].semilogy(x, z['rho'], ls, color='k', label=lab); A[1].semilogy(x, z['T'], ls, color='k', label=lab)
    A[2].semilogy(x, 1e-5 * z['vr_rms_fl'], ls, color='C0', label='v_r\' ' + lab); A[2].semilogy(x, 1e-5 * z['vt_rms'], ls, color='C1', label='v_th ' + lab); A[2].semilogy(x, 1e-5 * z['vp_rms'], ls, color='C2', label='v_ph ' + lab)
    A[3].semilogy(x, z['mach_t_rms'], ls, color='C0', label='M (gas+rad c_s) rms ' + lab); A[3].semilogy(x, z['mach_t_max'], ls, color='C3', label='max ' + lab)
    A[4].plot(x, S['nab'] - S['nad'], ls, color='k', label=lab); A[5].plot(x, S['N2'] / (S['g'] / S['Hp']), ls, color='k', label=lab)
    A[6].plot(x, z['Lrad_com'] / LS, ls, color='C0', label='comoving L_rad ' + lab); A[6].plot(x, z['Lenth'] / LS, ls, color='C1', label='L_enth (gas+rad) ' + lab); A[6].plot(x, z['Lkin'] / LS, ls, color='C2', label='L_KE ' + lab)
    A[7].semilogy(x, z['kapR'], ls, color='k', label='kappa_R ' + lab)
for a, ti, yl in zip(A, ('density', 'temperature', 'rms velocity (v_r about shell mean)', 'Mach number', 'superadiabaticity', 'N^2 H_P/g (shell means)', 'luminosities / 6.618e38', 'shell-mean kappa_R'), ('rho [g/cm3]', 'T [K]', 'km/s', 'M', 'nabla - nabla_ad', 'N2/(g/H)', 'L/L*', 'cm2/g')):
    a.set_title(ti); a.set_xlabel('r [R_sun]'); a.set_ylabel(yl); a.legend(fontsize=7); a.grid(True, lw=0.4); marks(a)
A[3].set_ylim(1e-5, 3)
A[4].axhline(0, color='grey', lw=0.6); A[5].axhline(0, color='grey', lw=0.6)
A[5].set_ylim(-1, 1); A[5].set_title('N^2 H_P/g (clipped to +-1)')
A[6].clear(); A[6].set_title('|L_enth|, |L_KE| / L* at tau > 1 (log; dotted = negative)'); A[6].set_xlabel('r [R_sun]')
for z, c0, lab in ((Z[iE], 0.4, 'early'), (zl, 1.0, 'last')):
    w = z['tau_sh'] > 1
    for k, c in (('Lenth', 'C1'), ('Lkin', 'C2')):
        q = np.where(w, z[k], np.nan) / LS; A[6].semilogy(x, np.where(q > 0, q, np.nan), '-', color=c, alpha=c0, label=k + ' ' + lab); A[6].semilogy(x, np.where(q < 0, -q, np.nan), ':', color=c, alpha=c0)
A[6].set_ylim(1e-9, 1e-2); A[6].set_xlim(7.5, 7.95); A[6].legend(fontsize=7); A[6].grid(True, lw=0.4); marks(A[6])
fig.suptitle('M35ZAMS (T2, rc4v): shell statistics; blue dotted = tau 1, red dashed = iron-peak shell (T %.2e K)' % zl['T'][iFe])
fig.savefig(os.path.join(OUT, 'fig1_profiles.png'), dpi=85); plt.close(fig)
f2, b2 = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
for z, ls, lab in ((Z[iE], '--', 'early'), (zl, '-', 'last')):
    for k, c in (('vr_rms_fl', 'C0'), ('vt_rms_fl', 'C1'), ('vp_rms_fl', 'C2')):
        b2[0].semilogy(x, 1e-5 * z[k + '_m1'], ls, color=c, label=k + ' ' + lab); b2[0].semilogy(x, 1e-5 * z[k], ls, color=c, alpha=0.3)
    b2[1].semilogy(x, z['mach_t_rms_m1'], ls, color='C0', label='rms ' + lab); b2[1].semilogy(x, z['mach_t_max_m1'], ls, color='C3', label='max ' + lab); b2[1].semilogy(x, z['mach_t_rms'], ls, color='C0', alpha=0.3); b2[1].semilogy(x, z['mach_t_max'], ls, color='C3', alpha=0.3)
    b2[2].semilogy(x, z['dT_rms_m1'], ls, color='k', label='rms T prime / T, ' + lab)
for a, ti in zip(b2, ('rms velocities, shell mean removed [km/s]', 'Mach (gas+rad c_s)', 'rms T fluctuation')):
    a.set_title(ti + '; cells tau >= 1 (faint: all cells)', fontsize=9); a.set_xlabel('r [R_sun]'); a.set_xlim(7.3, 8.05); a.legend(fontsize=7); a.grid(True, lw=0.4); marks(a)
b2[0].set_ylim(1e-6, 300); b2[1].set_ylim(1e-6, 5)
f2.savefig(os.path.join(OUT, 'fig4_tau1_excluded.png'), dpi=85); plt.close(f2)
fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
for j, (a, b) in enumerate(bands):
    ax[0].semilogy(t / 1e3, KF[:, j], label='%.3f-%.3f' % (x[a], x[min(b, 383)]), lw=1)
ax[0].set_title('KE_fl by radial shell'); ax[0].legend(fontsize=6, ncol=2)
ax[1].semilogy(t / 1e3, KTOT, 'k-', label='fluctuation total (v_r\' + lateral)'); ax[1].semilogy(t / 1e3, KLM, 'r--', label='shell-mean lateral flow'); ax[1].semilogy(t / 1e3, KRF, 'b:', label='radial fluctuation only'); ax[1].semilogy(t / 1e3, KMT, 'k--', label='mean radial flow total'); ax[1].set_title('total KE (4 pi scaled)'); ax[1].legend(fontsize=8)
ax[2].semilogy(h[:, 0] / 1e3, h[:, 7], label='1-KE (hst)'); ax[2].semilogy(h[:, 0] / 1e3, h[:, 8], label='2-KE'); ax[2].semilogy(h[:, 0] / 1e3, h[:, 9], label='3-KE'); ax[2].set_title('hydro.hst wedge KE'); ax[2].legend(fontsize=8)
for a in ax:
    a.set_xlabel('t [ks]'); a.grid(True, lw=0.4)
fig.savefig(os.path.join(OUT, 'fig2_ke_time.png'), dpi=90); plt.close(fig)
# slices at the iron peak, two times
import red as RR
fig, ax = plt.subplots(2, 4, figsize=(19, 8.4), constrained_layout=True)
for q, (n, lab) in enumerate(((int(Z[iE]['n']), 't=%.0f s' % t[iE]), (int(zl['n']), 't=%.0f s' % t[iL]))):
    hh = R.bc.read_binary(R.D + 'hydro_w.%05d.bin' % n)
    x1f, x1v, ths, phs, W, asm = R.grid(hh)
    rho = asm('dens', hh); vr = asm('velx', hh); ei = asm('eint', hh); vv = asm('vely', hh); wv = asm('velz', hh)
    T = (2 / 3.) * ei / rho * R.TU
    i = iFe; vb = vr[:, :, i] - np.average(vr[:, :, i], weights=W * rho[:, :, i]); Tb = T[:, :, i] / np.average(T[:, :, i], weights=W) - 1
    cc = np.average(vb * Tb, weights=W) / np.sqrt(np.average(vb ** 2, weights=W) * np.average(Tb ** 2, weights=W)); dens = rho[:, :, i] / np.average(rho[:, :, i], weights=W) - 1
    cr = np.average(vb * dens, weights=W) / np.sqrt(np.average(vb ** 2, weights=W) * np.average(dens ** 2, weights=W))
    sk = np.average(vb ** 3, weights=W) / np.average(vb ** 2, weights=W) ** 1.5; fd = float(W[vb < 0].sum() / W.sum())
    p('%s iron-peak shell: corr(v_r\', T\') = %.3f, corr(v_r\', rho\') = %.3f, skewness of v_r\' = %.2f, area fraction of downflow = %.2f, rms T\'/T = %.2e, rms rho\'/rho = %.2e, v_r\' min/max = %.2f/%.2f km/s' % (lab, cc, cr, sk, fd, np.sqrt(np.average(Tb ** 2, weights=W)), np.sqrt(np.average(dens ** 2, weights=W)), 1e-5 * vb.min(), 1e-5 * vb.max()))
    d_ = vr[:, :, i] - np.average(vr[:, :, i], weights=W); T_ = T[:, :, i] / np.average(T[:, :, i], weights=W) - 1
    vth_ = vv[:, :, i] / 1e5; s = max(1e-9, np.percentile(abs(d_), 99)) / 1e5
    xx = np.degrees(np.array(phs)); yy = np.degrees(np.array(ths))
    for c, (f_, ti, vm, cmap) in enumerate(((vr[:, :, i].T * 1e-5 - np.average(vr[:, :, i], weights=W) * 1e-5, "v_r' at iron peak [km/s]", s, 'RdBu_r'), (T_.T * 100, "T'/<T> [%]", max(1e-9, np.percentile(abs(T_), 99) * 100), 'PuOr_r'), (vth_.T, 'v_theta [km/s]', max(1e-9, np.percentile(abs(vth_), 99)), 'RdBu_r'), (np.log10(rho[:, :, i]).T, 'log rho', None, 'viridis'))):
        a = ax[q, c]; kw = dict(cmap=cmap, shading='auto'); kw.update(dict(vmin=-vm, vmax=vm) if vm else {})
        im = a.pcolormesh(xx, yy, f_, **kw); fig.colorbar(im, ax=a); a.set_title('%s, %s (r=%.3f R_sun)' % (ti, lab, x[i]), fontsize=9); a.set_xlabel('phi [deg]'); a.set_ylabel('theta [deg]')
fig.suptitle('M35ZAMS lateral slices at the iron-peak shell (T ~ %.2e K); colour limits = 99%% of |.| per panel' % zl['T'][iFe])
fig.savefig(os.path.join(OUT, 'fig3_slices_iron.png'), dpi=85); plt.close(fig)
open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(TXT) + '\n')
p('figs done')
