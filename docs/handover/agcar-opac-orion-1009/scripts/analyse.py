"""Assemble the production tables, validate against Ferguson 2005 / TOPS, compare with
viper's ext2, frequency-grid convergence and sensitivity; plots (cgs).
usage: python analyse.py  -> docs/handover/agcar-opac-orion-1009/{tables,plots}/"""
import glob, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/orion/ptmp/jinma/agcar_opac_1009/scripts')
import refs

W = '/orion/ptmp/jinma/agcar_opac_1009/work/out/'
O = '/orion/ptmp/jinma/agcar_opac_1009/repo/docs/handover/agcar-opac-orion-1009/'
EXT2 = ('/orion/ptmp/jinma/agcar_opac_1009/ext2/docs/handover/agcar-opac-1009/'
        '%s_ext2_gs98_x0.36_z0.02.txt')
os.makedirs(O + 'tables', exist_ok=True)
os.makedirs(O + 'plots', exist_ok=True)
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.color': '#e4e4e0',
                     'grid.linewidth': 0.5, 'lines.linewidth': 1.6})
REP = []


def say(s=''):
    print(s)
    REP.append(s)


# ---------------------------------------------------------------- production tables
LT = np.round(np.arange(3.45, 4.3001, 0.025), 3)
d0 = np.load(W + 'T%.3f_R1e+06_prod.npz' % LT[0])
LD = d0['lr']
KEYS = ['kP', 'kP_num', 'kR', 'kR_cont', 'kes', 'ne', 'kP_lines', 'kP_bf', 'kP_ff', 'kP_hm',
        'kP_mol']
P = {k: np.zeros((len(LT), len(LD))) for k in KEYS}
for i, t in enumerate(LT):
    d = np.load(W + 'T%.3f_R1e+06_prod.npz' % t)
    for k in KEYS:
        P[k][i] = d[k]
lP, lR = np.log10(P['kP']), np.log10(P['kR'])
say('# numerical (R=1e6 opacity sampling) vs exact line-sum Planck: max |dlog| = %.4f'
    % np.abs(np.log10(P['kP_num']) - lP).max())


def write_table(fn, K, what):
    with open(fn, 'w') as f:
        f.write(f'# AthenaK opacity table, log10 kappa_{what} [cm^2/g] on (log10 T[K], '
                'log10 rho[g/cm3])\n')
        f.write('# source: orion 10-09 LTE calculation (scripts/run_T.py): Kurucz atomic '
                'lines (CD-ROM 1 for the Fe group, gfall08oct17 otherwise), Verner+96 / '
                'hydrogenic bf, ff, H-, CO/CN/OH/SiO/H2O/TiO (DACE+FastChem)\n')
        f.write('# GS98 metals, X=0.360 Y=0.620 Z=0.020; '
                + ('ABSORPTION ONLY (no scattering)\n' if what == 'P' else
                   'incl. electron (Thomson) + Rayleigh scattering\n'))
        f.write('# VALID BOX = the whole table: log T 3.45..4.30, log rho -21..-12\n')
        f.write('# grid: nT nD lTmin dlT lDmin dlD\n')
        f.write('# %d %d %.3f 0.025 %.2f 0.05\n' % (len(LT), len(LD), LT[0], LD[0]))
        f.write('# then nT*nD rows, T slowest: log10 kappa [cm^2/g]\n')
        for v in K.ravel():
            f.write('%.5f\n' % v)


write_table(O + 'tables/planck_lowrho_orion_gs98_x0.36_z0.02.txt', lP, 'P')
write_table(O + 'tables/rosseland_lowrho_orion_gs98_x0.36_z0.02.txt', lR, 'R')
# component summary on a coarse subset
with open(O + 'tables/planck_components.txt', 'w') as f:
    f.write('# log T, log rho, log kP, fractions of kP: lines bf ff H- molecules; log kR, '
            'log kes, log kR_cont(no lines), log n_e\n')
    for i in range(0, len(LT), 2):
        for j in range(0, len(LD), 10):
            k = P['kP'][i, j]
            f.write('%.3f %6.2f %7.3f  %.3f %.2e %.2e %.2e %.2e  %7.3f %7.3f %7.3f %6.2f\n' % (
                LT[i], LD[j], lP[i, j], P['kP_lines'][i, j]/k, P['kP_bf'][i, j]/k,
                P['kP_ff'][i, j]/k, P['kP_hm'][i, j]/k, P['kP_mol'][i, j]/k, lR[i, j],
                np.log10(P['kes'][i, j]), np.log10(P['kR_cont'][i, j]),
                np.log10(P['ne'][i, j])))
fm = P['kP_mol']/P['kP']
im_, jm_ = np.unravel_index(fm.argmax(), fm.shape)
say('# molecules: max fraction of kP %.3g at (log T, log rho) = (%.3f, %.2f)'
    % (fm.max(), LT[im_], LD[jm_]))
for thr in [0.01, 0.1]:
    m = (P['kP_mol']/P['kP']) > thr
    if m.any():
        ii, jj = np.where(m)
        say('#   molecules > %g of kP at log T <= %.3f, log rho >= %.2f (%d cells)'
            % (thr, LT[ii].max(), LD[jj].min(), m.sum()))
mh = (P['kP_hm'] + P['kP_bf'] + P['kP_ff'])/P['kP']
say('# continuum (bf+ff+H-) max fraction of kP: %.3g' % mh.max())
say('# Rosseland: max kR/kR_cont (line contribution) %.3f; kR/kes range %.2f..%.2f'
    % ((P['kR']/P['kR_cont']).max(), (P['kR']/P['kes']).min(), (P['kR']/P['kes']).max()))

# ---------------------------------------------------------------- validation
TT, DD = np.meshgrid(LT, LD, indexing='ij')
v2 = '/orion/ptmp/jinma/agcar_opac_1009/ext2/docs/handover/agcar-opac-1009/sources/'
for kind, mine in [('planck', lP), ('rosseland', lR)]:
    F = refs.Ferg(kind)(TT, DD)
    Tp = refs.Tops(kind)(TT, DD)
    for nm, ref in [('Ferguson05', F), ('TOPS', Tp)]:
        m = np.isfinite(ref)
        dd = mine[m] - ref[m]
        say(f'# {kind} vs {nm}: N={m.sum()} median {np.median(dd):+.3f} dex, '
            f'p10 {np.percentile(dd, 10):+.3f} p90 {np.percentile(dd, 90):+.3f}, '
            f'max |d| {np.abs(dd).max():.3f}; |d|<0.2 dex in {np.mean(abs(dd) < 0.2)*100:.0f} %')
        # per-T row summary
        rows = []
        for i in range(0, len(LT), 2):
            mm = m[i]
            if mm.sum() == 0:
                continue
            x = mine[i, mm] - ref[i, mm]
            rows.append('%.3f: n %2d med %+.2f min %+.2f max %+.2f (log rho %.2f..%.2f)'
                        % (LT[i], mm.sum(), np.median(x), x.min(), x.max(), LD[mm].min(),
                           LD[mm].max()))
        with open(O + f'tables/validation_{kind}_{nm}.txt', 'w') as f:
            f.write(f'# log10 kappa_{kind[0].upper()}(this work) - {nm}, per log T row over '
                    'the overlap\n' + '\n'.join(rows) + '\n')
    # Ferguson vs TOPS in the same overlap (the reference disagreement)
    m = np.isfinite(F) & np.isfinite(Tp)
    if m.any():
        dd = F[m] - Tp[m]
        say(f'# {kind}: Ferguson05 - TOPS on the common overlap: median {np.median(dd):+.3f},'
            f' p10 {np.percentile(dd, 10):+.3f} p90 {np.percentile(dd, 90):+.3f}')

# ---------------------------------------------------------------- ext2
E = {}
for kind in ['planck', 'rosseland']:
    t, dgrid, K = refs.read_repo(EXT2 % kind)
    it = [int(np.argmin(abs(t - x))) for x in LT]
    jd = [int(np.argmin(abs(dgrid - x))) for x in LD]
    assert np.allclose(t[it], LT) and np.allclose(dgrid[jd], LD)
    E[kind] = K[np.ix_(it, jd)]
with open(O + 'tables/ext2_compare.txt', 'w') as f:
    f.write('# log10 kappa(this work) - log10 kappa(ext2) on the ext2 nodes, per log T row; '
            'columns: log rho = -21 -20 -19 -18 -17 -16 -15 -14 -13 -12\n')
    for kind, mine in [('planck', lP), ('rosseland', lR)]:
        dd = mine - E[kind]
        f.write(f'## {kind}\n')
        for i in range(len(LT)):
            f.write('%.3f ' % LT[i] + ' '.join('%+6.2f' % dd[i, j] for j in range(0, len(LD), 20))
                    + '\n')
        for lo, hi in [(-21, -15), (-15, -12)]:
            m = (DD >= lo) & (DD <= hi)
            say(f'# {kind} this - ext2, log rho {lo}..{hi}: median {np.median(dd[m]):+.2f}, '
                f'p10 {np.percentile(dd[m], 10):+.2f}, p90 {np.percentile(dd[m], 90):+.2f},'
                f' max |d| {np.abs(dd[m]).max():.2f}')
        m = (DD <= -15) & (TT >= 3.5)
        say(f'# {kind} this - ext2 at rho 1e-20..1e-15 (log T >= 3.5): '
            f'min {dd[(DD >= -20) & (DD <= -15)].min():+.2f} max {dd[(DD >= -20) & (DD <= -15)].max():+.2f}')

# ---------------------------------------------------------------- convergence, sensitivity
LTc = np.round(np.arange(3.45, 4.3001, 0.05), 2)


def coarse(tag, R='1e+05'):
    out = {}
    for k in ['kP', 'kP_num', 'kR']:
        out[k] = np.array([np.load(W + 'T%.3f_R%s%s.npz' % (t, R, tag))[k] for t in LTc])
    return out


ref6 = {k: v[np.ix_([int(np.argmin(abs(LT - x))) for x in LTc], range(0, len(LD), 5))]
        for k, v in P.items()}
with open(O + 'tables/convergence.txt', 'w') as f:
    f.write('# frequency-grid convergence: opacity sampling on a log-uniform grid 300..4e5 '
            'cm^-1 (33 um..25 nm) with R = nu/dnu; columns: R, max and median |dlog| vs '
            'R=1e6 of kP_num (sampled), kR; exact line-sum kP is grid independent\n')
    for R in ['1e+04', '3e+04', '1e+05', '3e+05']:
        c = coarse('', R)
        a = np.abs(np.log10(c['kP_num']/ref6['kP_num']))
        b = np.abs(np.log10(c['kR']/ref6['kR']))
        e = np.abs(np.log10(c['kP']/ref6['kP']))
        line = (f'R {R}: kP_num max {a.max():.4f} med {np.median(a):.4f} | kR max {b.max():.4f}'
                f' med {np.median(b):.4f} | kP exact max {e.max():.2e}')
        f.write(line + '\n')
        say('# conv ' + line)
base = coarse('')
with open(O + 'tables/sensitivity.txt', 'w') as f:
    f.write('# dlog10 vs the R=1e5 baseline (mrg lines, WCAP 0.01, xi 2 km/s, NMAX_H 30, '
            'EPS 1e-10): max |d| and median d over the coarse grid\n')
    for tag, lab in [('_gf08', 'gfall08oct17 only (no predicted Fe-group lines)'),
                     ('_wcap003', 'line wings cut at 0.3 % of nu'),
                     ('_wcap03', 'line wings cut at 3 % of nu'),
                     ('_xi0', 'no microturbulence'), ('_nh15', 'H I/He II n_max 15'),
                     ('_nh60', 'H I/He II n_max 60'), ('_eps12', 'line skip threshold 1e-12')]:
        c = coarse(tag)
        for k in ['kP', 'kR']:
            dd = np.log10(c[k]/base[k])
            i, j = np.unravel_index(np.abs(dd).argmax(), dd.shape)
            line = (f'{lab:48s} {k}: max |d| {np.abs(dd).max():.3f} at (log T {LTc[i]:.2f}, '
                    f'log rho {-21 + 0.25*j:.2f}), median {np.median(dd):+.4f}')
            f.write(line + '\n')
            say('# sens ' + line)

# ---------------------------------------------------------------- plots
ext = [LD[0], LD[-1], LT[0], LT[-1]]
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
for a, Z, ttl, cm, lim in [(ax[0], lP, r'this work: log $\kappa_P$ [cm$^2$ g$^{-1}$]', 'viridis',
                            None),
                           (ax[1], E['planck'], r'ext2: log $\kappa_P$', 'viridis', None),
                           (ax[2], lP - E['planck'], r'this work $-$ ext2 [dex]', 'RdBu_r',
                            (-2, 2))]:
    if lim is None:
        lim = (np.floor(min(lP.min(), E['planck'].min())), np.ceil(max(lP.max(), E['planck'].max())))
    im = a.imshow(Z, origin='lower', aspect='auto', extent=ext, cmap=cm, vmin=lim[0], vmax=lim[1])
    plt.colorbar(im, ax=a)
    a.set_title(ttl)
    a.set_xlabel(r'log $\rho$ [g cm$^{-3}$]')
    a.grid(False)
ax[0].set_ylabel('log T [K]')
for a in ax:
    a.plot(3*LT - 26, LT, color='w', lw=1, ls='--')      # Ferguson log R = -8 edge
    a.axvline(-14, color='w', lw=1, ls=':')
fig.suptitle(r'LTE Planck mean, AG Car X 0.36 Z 0.02 (dashed: Ferguson log R=$-8$ edge; '
             r'dotted: TOPS $\rho$ = 1e-14)')
fig.savefig(O + 'plots/planck_map_vs_ext2.png', dpi=130)
plt.close(fig)

fig, ax = plt.subplots(2, 3, figsize=(12, 6.5), constrained_layout=True, sharex=True)
for a, t in zip(ax.ravel(), [3.45, 3.6, 3.8, 3.95, 4.1, 4.3]):
    i = int(np.argmin(abs(LT - t)))
    rho = 10**LD
    a.semilogx(rho, lP[i], color=C[0], label='this work (LTE)')
    a.semilogx(rho, E['planck'][i], color=C[3], ls='--', label='ext2 (viper)')
    a.semilogx(rho, refs.Tops('planck')(LT[i], LD), color=C[1], marker='o', ms=3, ls='none',
               label='TOPS')
    a.semilogx(rho, refs.Ferg('planck')(LT[i], LD), color=C[2], marker='s', ms=3, ls='none',
               label='Ferguson 2005')
    a.set_title('T = %.0f K' % 10**LT[i])
    a.set_ylabel(r'log $\kappa_P$ [cm$^2$ g$^{-1}$]')
for a in ax[1]:
    a.set_xlabel(r'$\rho$ [g cm$^{-3}$]')
ax[0, 0].legend(frameon=False, fontsize=8)
fig.savefig(O + 'plots/planck_cuts.png', dpi=130)
plt.close(fig)

fig, ax = plt.subplots(2, 3, figsize=(12, 6.5), constrained_layout=True, sharex=True)
for a, t in zip(ax.ravel(), [3.45, 3.6, 3.8, 3.95, 4.1, 4.3]):
    i = int(np.argmin(abs(LT - t)))
    rho = 10**LD
    a.semilogx(rho, lR[i], color=C[0], label='this work (LTE)')
    a.semilogx(rho, np.log10(P['kes'][i]), color='#8a8a85', lw=1, ls=':', label='Thomson only')
    a.semilogx(rho, E['rosseland'][i], color=C[3], ls='--', label='ext2 (viper)')
    a.semilogx(rho, refs.Tops('rosseland')(LT[i], LD), color=C[1], marker='o', ms=3, ls='none',
               label='TOPS')
    a.semilogx(rho, refs.Ferg('rosseland')(LT[i], LD), color=C[2], marker='s', ms=3, ls='none',
               label='Ferguson 2005')
    a.set_title('T = %.0f K' % 10**LT[i])
    a.set_ylabel(r'log $\kappa_R$ [cm$^2$ g$^{-1}$]')
for a in ax[1]:
    a.set_xlabel(r'$\rho$ [g cm$^{-3}$]')
ax[0, 0].legend(frameon=False, fontsize=8)
fig.savefig(O + 'plots/rosseland_cuts.png', dpi=130)
plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
Rs = [1e4, 3e4, 1e5, 3e5]
for k, a, lab in [('kP_num', ax[0], r'sampled $\kappa_P$'), ('kR', ax[1], r'$\kappa_R$')]:
    mx, md = [], []
    for R in ['1e+04', '3e+04', '1e+05', '3e+05']:
        c = coarse('', R)
        x = np.abs(np.log10(c[k]/ref6[k]))
        mx.append(x.max()), md.append(np.median(x))
    a.loglog(Rs, mx, color=C[0], marker='o', label='max over grid')
    a.loglog(Rs, md, color=C[1], marker='s', label='median')
    a.set_xlabel(r'$R=\nu/\Delta\nu$')
    a.set_ylabel(r'|$\Delta$ log| vs R = 1e6 [dex]')
    a.set_title(lab)
ax[0].legend(frameon=False)
fig.savefig(O + 'plots/convergence.png', dpi=130)
plt.close(fig)

fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
cs = ax.contourf(LD, LT, np.log10(np.maximum(P['kP_mol']/P['kP'], 1e-6)),
                 levels=np.arange(-6, 0.01, 1), cmap='Blues')
plt.colorbar(cs, ax=ax, label=r'log (molecular share of $\kappa_P$)')
ax.set_xlabel(r'log $\rho$ [g cm$^{-3}$]')
ax.set_ylabel('log T [K]')
ax.set_ylim(3.45, 3.7)
ax.set_title('where molecules (CO, CN, OH, SiO, H2O, TiO) matter')
fig.savefig(O + 'plots/molecule_share.png', dpi=130)
plt.close(fig)
open(O + 'tables/summary.txt', 'w').write('\n'.join(REP) + '\n')
