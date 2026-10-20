#!/usr/bin/env python3
"""ACC3D page 10-20: TRACER section (READ ONLY on the run). usage: mk_tracer.py <imgdir>  -> <imgdir>/trc.json + figures.
Tracers (deck v11amwtr, <hydro>/nscalars = 2): s0 = stream gas (set to 1 by the nozzle reset), s1 = gas that entered through
the other user boundaries. ryper.trc.hst: Ms0 Ms1 (domain integrals of rho s), cumulative Nz0 (s0 written by the nozzle), Nz1,
Sk0 Sk1 (removed by the sinks), Fi0 Fi1 (out through r_in), Fo0 Fo1 (out through r_out + theta walls), Bin (mass in through
the non-nozzle boundaries); budget Ms0 - Ms0(0) = Nz0 - Sk0 - Fi0 - Fo0.
Per selected stream dump (NTRC times, cached in $OUT/data/trc/): phi-summed s0 maps, the s0 mass profile in r, and the
s0-weighted gas entropy / optical depth of the accreted gas (eload.py: true faces, EOS-table T and s). Without trc.hst or
without s_00 in the dumps the section is skipped (trc.json says why)."""
import os
import sys
import json
import numpy as np
import common as C

IMG = sys.argv[1]
CD = os.path.join(C.OUT, 'data', 'trc'); os.makedirs(CD, exist_ok=True)
NTRC = int(C.E('NTRC', '4'))
R_TOP = 10.3          # top of the region counted as "on the star" (above: stream + atmosphere), as entropy.py (run11)
SB = np.arange(30, 66, 0.25)
TB = np.arange(-3, 6.01, 0.1)    # log10 tau bins
res = dict(present=False)


def save():
    json.dump(res, open(os.path.join(IMG, 'trc.json'), 'w'), indent=1)


T = C.hstcat('trc')
if T is None:
    res['why'] = 'no ryper.trc.hst in %s (run without tracers)' % C.stream_dirs_label()
    save(); print(res['why']); sys.exit(0)
plt = C.mpl()
t = T['time']
res['present'] = True
on = T['Nz0'] > 0
t_on = float(t[np.argmax(on)]) if on.any() else C.TS
fill = float(T['Nz0'][np.argmax(on)]) if on.any() else 0.0   # the nozzle cells filled at switch-on (first row)
bud0 = (T['Ms0'] - T['Ms0'][0]) - (T['Nz0'] - T['Sk0'] - T['Fi0'] - T['Fo0'])
bud1 = (T['Ms1'] - T['Ms1'][0]) - (-T['Nz1'] - T['Sk1'] - T['Fi1'] - T['Fo1'])
late = t > max(t[-1] - 0.05, t_on + 0.005)
if late.sum() > 2:
    p = np.polyfit(t[late], T['Nz0'][late], 1); q = np.polyfit(t[late], T['Ms0'][late], 1)
    res['Nz0_rate'] = float(p[0]); res['Ms0_rate'] = float(q[0])
res.update(t=float(t[-1]), t_on=t_on, fill=fill, Ms0=float(T['Ms0'][-1]), Nz0=float(T['Nz0'][-1]), Sk0=float(T['Sk0'][-1]),
           Fi0=float(T['Fi0'][-1]), Fo0=float(T['Fo0'][-1]), Ms1=float(T['Ms1'][-1]), Bin=float(T['Bin'][-1]),
           Fi1=float(T['Fi1'][-1]), Fo1=float(T['Fo1'][-1]), Sk1=float(T['Sk1'][-1]), Nz1=float(T['Nz1'][-1]),
           nominal=float(fill + C.MDOT * max(t[-1] - t_on, 0)), bud0=float(np.max(np.abs(bud0))), bud1=float(np.max(np.abs(bud1))),
           s0max=float(np.max(T['s0max'])), s0min=float(np.min(T['s0min'])), smax=float(np.max(T['smax'])))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
a = ax[0]
a.plot(t, T['Ms0'], color=C.C1, lw=1.8, label='Ms0: s0 mass in the domain')
a.plot(t, T['Nz0'], color=C.C2, lw=1, label='Nz0: written by the nozzle (cum.)')
a.plot(t[t >= t_on], fill + C.MDOT * (t[t >= t_on] - t_on), color='#a3adb8', ls='--', lw=.8, label='nozzle fill + %.2f x (t - t$_{on}$)' % C.MDOT)
a.plot(t, T['Fi0'], color=C.C3, lw=1, label='Fi0: out through r$_{in}$ (cum.)')
a.plot(t, T['Fo0'], color=C.C4, lw=1, label=r'Fo0: out through r$_{out}$ + $\theta$ walls')
a.plot(t, T['Sk0'], color=C.C5, lw=1, label='Sk0: removed by the sinks')
a.plot(t, T['Ms1'], color='#e8ecf1', lw=1, ls=':', label='Ms1: boundary-inflow tracer mass')
a.set_title('tracer masses [code mass]'); a.legend(fontsize=7)
a = ax[1]
tm, r0 = C.rate(t, T['Nz0'], 9); _, r1 = C.rate(t, T['Ms0'], 9)
a.plot(tm, r0, color=C.C2, lw=1, label='dNz0/dt (nozzle)'); a.plot(tm, r1, color=C.C1, lw=1, label='dMs0/dt (domain)')
a.axhline(C.MDOT, color='#a3adb8', ls='--', lw=.8, label='nominal %.2f' % C.MDOT)
a.set_xlim(t_on - 0.002, None); a.set_ylim(-50, 1.6 * C.MDOT)
a.set_title('s0 rates [code mass / code time] (9-row mean)'); a.legend(fontsize=7)
a = ax[2]
a.plot(t, bud0, color=C.C1, label='s0: dMs0 - (Nz0 - Sk0 - Fi0 - Fo0)')
a.plot(t, bud1, color=C.C2, label='s1: dMs1 + (Nz1 + Sk1 + Fi1 + Fo1)')
a.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3)); a.set_title('tracer budget residuals [code mass]'); a.legend(fontsize=7)
for a in ax:
    a.axvline(C.TS, color='#a3adb8', lw=.6, ls=':'); a.grid(alpha=.5); a.set_xlabel('t [code]')
fig.suptitle(f'Tracer history (ryper.trc.hst), {C.RUN_LABEL} ({C.SITE})', color='#e8ecf1')
fig.tight_layout(); fig.savefig(IMG + '/fig_trc_hst.png', dpi=80); plt.close(fig)
print('trc hst done', flush=True)

# ---------------- per-dump analysis ----------------
L = C.dumps('hydro_w')
sm = [x for x in L if x[2] == 'stream' and x[0] > t_on + float(C.E('TRC_TMIN', '0.005'))]
stt = [x for x in L if x[2] == 'static']
if not sm:
    res['dumps'] = []; res['why_dumps'] = 'no stream-phase dump yet'
    save(); print(res['why_dumps']); sys.exit(0)
tg = np.linspace(sm[0][0], sm[-1][0], NTRC) if len(sm) > 1 else [sm[-1][0]]
sel = sorted({min(range(len(sm)), key=lambda i: abs(sm[i][0] - x)) for x in tg} | {len(sm) - 1})
sel = [sm[i] for i in sel]


def wpct(x, w, q):
    if w.sum() <= 0:
        return [float('nan')] * len(q)
    o = np.argsort(x); c = np.cumsum(w[o]); c = c / c[-1]
    return [float(v) for v in np.interp(np.asarray(q) / 100.0, c, x[o])]


def analyse(f, ref=False):
    key = os.path.join(CD, '%s_%05d.npz' % ('ref' if ref else 'trc', C.dnum(f)))
    if os.path.exists(key) and os.path.getmtime(key) > os.path.getmtime(f) and not C.FORCE:
        z = np.load(key, allow_pickle=True); return {k: (z[k].item() if z[k].dtype == object else z[k]) for k in z.files}
    import eload
    D = eload.load(f, tau=True)
    R3 = np.broadcast_to(D['r'][None, None, :], D['rho'].shape)
    dm = D['dm']
    o = dict(t=D['t'], r=D['r'], rf=D['rf'], thf=D['thf'], th=D['th'], ph=D['ph'])
    o['M_tr'] = dm.sum(axis=0)
    o['s_shell'] = (dm * D['s']).sum(axis=(0, 1)) / dm.sum(axis=(0, 1))
    lat = np.abs(90 - np.degrees(D['th']))
    eqm = (lat < 10)[None, :, None]
    o['s_shell_eq'] = (dm * D['s'] * eqm).sum(axis=(0, 1)) / (dm * eqm).sum(axis=(0, 1))
    if 'tau' in D:
        o['tau_shell'] = np.exp((dm * np.log(np.maximum(D['tau'], 1e-30))).sum(axis=(0, 1)) / dm.sum(axis=(0, 1)))
    surf = (R3 > 8.8) & (R3 < 9.5)
    o['H_surf'] = np.histogram(D['s'][surf], SB, weights=dm[surf])[0]
    if 's_00' in D:
        s0 = np.clip(D['s_00'], 0, 1); w = dm * s0
        o['Ms0_tr'] = w.sum(axis=0)
        o['Ms0_r'] = w.sum(axis=(0, 1)); o['M_r'] = dm.sum(axis=(0, 1))
        rr = (R3 > 9.3) & (R3 < 10.5)
        kc = int(np.argmax((w * rr).sum(axis=(1, 2))))
        o['kc'] = kc; o['s0_cut'] = s0[kc].astype(np.float32)
        star = R3 < R_TOP
        ws = w[star]; x = D['s'][star]
        o['acc'] = dict(M=float(ws.sum()), s=float((ws * x).sum() / max(ws.sum(), 1e-300)),
                        srad=float((ws * D['srad'][star]).sum() / max(ws.sum(), 1e-300)),
                        sradE=float((ws * D['sradE'][star]).sum() / max(ws.sum(), 1e-300)) if 'sradE' in D else None,
                        pct=wpct(x, ws, [10, 50, 90]),
                        r=float((ws * R3[star]).sum() / max(ws.sum(), 1e-300)),
                        rpct=wpct(R3[star], ws, [10, 50, 90]),
                        tau=float(np.exp((ws * np.log(np.maximum(D['tau'][star], 1e-30))).sum() / max(ws.sum(), 1e-300))) if 'tau' in D else None,
                        taupct=wpct(np.log10(np.maximum(D['tau'][star], 1e-30)), ws, [10, 50, 90]) if 'tau' in D else None,
                        M_rmeas=float(w[R3 < C.RMEAS].sum()), M_racc=float(w[R3 < C.RACC].sum()), Mtot=float(w.sum()))
        o['H_acc'] = np.histogram(x, SB, weights=ws)[0]
        if 'tau' in D:
            o['Ht_acc'] = np.histogram(np.log10(np.maximum(D['tau'][star], 1e-30)), TB, weights=ws)[0]
        st_ = ~star
        o['stream'] = dict(M=float(w[st_].sum()), s=float((w[st_] * D['s'][st_]).sum() / max(w[st_].sum(), 1e-300)))
        if 's_01' in D:
            w1 = dm * np.clip(D['s_01'], 0, 1); o['Ms1_r'] = w1.sum(axis=(0, 1))
    np.savez(key, **{k: (np.array(v, dtype=object) if isinstance(v, dict) else v) for k, v in o.items()})
    return o


A = [analyse(f) for _, f, _ in sel]
if 'Ms0_r' not in A[-1]:
    res['dumps'] = []; res['why_dumps'] = 'the hydro_w dumps carry no s_00 (tracers off in the dump)'
    save(); print(res['why_dumps']); sys.exit(0)
R0 = analyse(stt[0][1], ref=True) if stt else None          # t = 0 star
RP = analyse(stt[-1][1], ref=True) if len(stt) > 1 else R0   # end of the static phase
res['dumps'] = [dict(t=float(a['t']), file=os.path.basename(f), acc=a['acc'], stream=a['stream'],
                     phi_impact=float(np.degrees(a['ph'][int(a['kc'])]))) for a, (_, f, _) in zip(A, sel)]
r = A[-1]['r']; rf = A[-1]['rf']; thf = A[-1]['thf']
cols = plt.cm.viridis(np.linspace(0.1, 0.95, len(A)))
lab = lambda a: 'stream +%.2f d (t %.3f)' % ((a['t'] - C.TS) * C.DAYS, a['t'])   # noqa: E731

# maps: phi-summed mass-weighted s0 fraction, and the cut through the impact azimuth (R_cyl - z, r <= d_L1)
ir = int(np.searchsorted(rf, C.DL1))
TT, RR = np.meshgrid(thf, rf[:ir + 1], indexing='ij')
X, Z = RR * np.sin(TT), RR * np.cos(TT)
zmax = C.DL1 * abs(np.cos(thf[0])) * 1.02
fig, ax = plt.subplots(len(A), 2, figsize=(12, 4.4 * len(A) + 0.6), squeeze=False)
for i, a in enumerate(A):
    fr = a['Ms0_tr'][:, :ir] / np.maximum(a['M_tr'][:, :ir], 1e-300)
    for j, (z, ttl) in enumerate(((fr, r'$\phi$-averaged, mass weighted $\langle s_0\rangle$'),
                                  (a['s0_cut'][:, :ir], r's$_0$ in the meridional cut through the impact ($\phi$ %.0f°)' % np.degrees(a['ph'][int(a['kc'])])))):
        b = ax[i, j]
        pc = b.pcolormesh(X, Z, np.log10(np.maximum(z, 1e-6)), cmap='viridis', vmin=-4, vmax=0, shading='flat', rasterized=True)
        b.add_patch(plt.Circle((0, 0), C.RACC, fill=False, ec='#e8ecf1', lw=.6, ls=':'))
        b.set_xlim(C.RIN, C.DL1); b.set_ylim(-zmax, zmax); b.set_aspect('equal')
        b.set_title(f'{lab(a)}: {ttl}', fontsize=9); b.set_ylabel('z [R$_\\odot$]')
        if i == len(A) - 1:
            b.set_xlabel(r'R$_{cyl}$ [R$_\odot$]  (dotted: R$_{acc}$ 9.0)')
        fig.colorbar(pc, ax=b, shrink=.9, pad=.01).set_label(r'log$_{10}$ s$_0$', fontsize=8)
fig.tight_layout(); fig.savefig(IMG + '/fig_trc_maps.png', dpi=72); plt.close(fig)

# profiles in r
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
dr = np.diff(rf)
for c, a in zip(cols, A):
    ax[0].semilogy(r, np.maximum(a['Ms0_r'] / dr, 1e-12), color=c, lw=1.4, label=lab(a))
    ax[1].semilogy(r, np.maximum(a['Ms0_r'] / np.maximum(a['M_r'], 1e-300), 1e-9), color=c, lw=1.4)
    cm_ = np.cumsum(a['Ms0_r']); ax[2].plot(r, cm_ / max(cm_[-1], 1e-300), color=c, lw=1.4)
ax[0].set_title(r'dM$_{s0}$/dr [code mass / R$_\odot$]'); ax[0].legend(fontsize=7)
ax[1].set_title('s0 mass fraction of the shell (all angles)'); ax[1].set_ylim(1e-8, 1.5)
ax[2].set_title(r'fraction of the domain s0 mass inside r')
for b in ax:
    b.set_xlim(8.0, C.DL1); b.grid(alpha=.5); b.set_xlabel(r'r [R$_\odot$]')
    b.axvline(C.RACC, color='0.6', ls='-.', lw=.8); b.axvline(R_TOP, color='0.6', ls=':', lw=.8)
fig.suptitle(f'Where the stream gas ends up ({C.RUN_LABEL}, {C.SITE}); dash-dot R$_{{acc}}$ 9.0, dotted r {R_TOP} (top of "on the star")', color='#e8ecf1')
fig.tight_layout(); fig.savefig(IMG + '/fig_trc_prof.png', dpi=80); plt.close(fig)

# entropy of the accreted (s0) gas
bc = 0.5 * (SB[1:] + SB[:-1])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
b = ax[0]


def hstep(b, h, **k):
    h = np.asarray(h, float); s = h.sum()
    if s > 0:
        b.step(bc, h / s / np.diff(SB), where='mid', **k)


if R0 is not None:
    hstep(b, R0['H_surf'], color='#e8ecf1', ls='--', lw=1, label='star t = 0, r 8.8-9.5 (all mass)')
for c, a in zip(cols, A):
    hstep(b, a['H_acc'], color=c, lw=1.4, label='s0 gas r < %.1f, %s' % (R_TOP, lab(a)))
b.set_xlim(32, 64); b.set_xlabel(r's$_{gas}$ [k$_B$/baryon]'); b.set_ylabel('mass fraction per unit s')
b.set_title('s0-weighted gas entropy of the accreted gas'); b.legend(fontsize=6.5)
b = ax[1]
if RP is not None:
    mm = (RP['r'] > 8.0) & (RP['r'] < 10.6)
    b.plot(RP['r'][mm], RP['s_shell_eq'][mm], color='#e8ecf1', lw=1.4, label='star, end of static phase (|lat| < 10°)')
for c, a in zip(cols, A):
    ac = a['acc']
    b.errorbar([ac['rpct'][1]], [ac['s']], xerr=[[ac['rpct'][1] - ac['rpct'][0]], [ac['rpct'][2] - ac['rpct'][1]]],
               yerr=[[ac['s'] - ac['pct'][0]], [ac['pct'][2] - ac['s']]], color=c, marker='o', ms=5, capsize=2, lw=1)
b.set_xlabel(r'r [R$_\odot$]'); b.set_ylabel(r's$_{gas}$ [k$_B$/baryon]'); b.set_ylim(33, 50)
b.set_title('accreted gas (dots: s0-weighted mean s, median r; bars p10-p90)', fontsize=9); b.legend(fontsize=7)
b = ax[2]
tc_ = 0.5 * (TB[1:] + TB[:-1])
for c, a in zip(cols, A):
    if 'Ht_acc' in a:
        h = a['Ht_acc']
        if h.sum() > 0:
            b.step(tc_, h / max(h.sum(), 1e-300), where='mid', color=c, lw=1.4, label=lab(a))
if RP is not None and 'tau_shell' in RP:
    for tv in (2 / 3., 100.):
        b.axvline(np.log10(tv), color='0.6', ls=':', lw=.8)
b.set_xlabel(r'log$_{10}\tau_R$ (radial, from r$_{out}$)'); b.set_ylabel('s0 mass fraction per bin')
b.set_title(r'optical depth of the accreted s0 gas (r < %.1f); dotted $\tau$ 2/3, 100' % R_TOP, fontsize=9); b.legend(fontsize=7)
for b in ax:
    b.grid(alpha=.5)
fig.tight_layout(); fig.savefig(IMG + '/fig_trc_entropy.png', dpi=80); plt.close(fig)
if RP is not None:
    # neutral-buoyancy radius: where the pre-stream star's equatorial shell entropy crosses the accreted value
    sa = A[-1]['acc']['s']; prof = RP['s_shell_eq']; ok = (RP['r'] > 8.0) & (RP['r'] < 9.7)
    ii = np.where(ok[:-1] & ok[1:] & (prof[:-1] < sa) & (prof[1:] >= sa))[0]
    if len(ii):
        i = ii[0]; fr = (sa - prof[i]) / (prof[i + 1] - prof[i])
        res['r_nb'] = float(RP['r'][i] + fr * (RP['r'][i + 1] - RP['r'][i]))
        if 'tau_shell' in RP:
            res['tau_nb'] = float(np.exp(np.log(RP['tau_shell'][i]) + fr * np.log(RP['tau_shell'][i + 1] / RP['tau_shell'][i])))
    for tv in (100.0, 2 / 3.):
        if 'tau_shell' in RP:
            k = int(np.argmin(np.abs(np.log(RP['tau_shell'] / tv))))
            res['star_s_tau%g' % round(tv, 2)] = [float(RP['r'][k]), float(RP['s_shell_eq'][k])]
save()
print(json.dumps(res, indent=1)[:3000])
print('tracer done')
