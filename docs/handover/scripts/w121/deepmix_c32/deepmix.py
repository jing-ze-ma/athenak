#!/usr/bin/env python3
"""deepmix.py: radial energy flux of the fluid per radial face of a cubed-sphere
deep_hot_jupiter_rt HYDRO run, split into mean-meridional / stationary-eddy / transient
parts, into enthalpy / kinetic / potential parts, into scales (coarse-graining), and
expressed as an effective entropy diffusivity K_eff compared with grid-scale estimates.

Built for the C32 WASP-121b productions (w121prod_0929/w1x, w10x); to be repeated
IDENTICALLY on C256 (question 09-30: is the upward deep enthalpy flux of 10x resolved
circulation or numerical mixing?).

INPUTS (all read in place, nothing is written outside --out)
  --rst   restart pattern, e.g. RUN/rst/dhj.%05d.rst (one file per restart; hydro only;
          the embedded input supplies grid, ap, grav, omega, rot_potential, etotgrav)
  --bin   ONE bin dump of the same run with the same MeshBlock layout (used only for the
          per-MeshBlock panel geometry: lat/lon of every column, via dhjcs.cs_geometry)
  --eost  EOS inversion table (dhjcs.EOS: eos_table.txt) of the run's composition
  --eosg  EOS grid (dclib.EOS: eos_grid.txt) of the run's composition (s, Gamma1)
  --win   windows "i0:i1[:stride],..." in restart INDEX (inclusive); rot = t/Prot from
          the restart header
  --prot  rotation period [s] (WASP-121b 1.101535e5)
  --out   output directory; writes deepmix_<label>.npz and deepmix_<label>.txt
  --ledger  optional ledger_*.txt (budget ledger.py output): its per-shell
          Fin(i)-Fout(i+1) column is cumulated from the wall (FXE(0)=Fbot=0) into the
          Riemann-flux energy flux per face, printed next to the resolved flux
  --nband zonal-mean latitude bands (default 30 = 6 deg; keep 30 on C256 so the split
          is defined on the SAME bands)
  --cg    horizontal coarse-graining factors (cells per side, must divide the MeshBlock
          nx2/nx3), default "1,2,4,8".  C32 cell = C256 x 8, so on C256 use
          "1,2,4,8,16,32,64" (MeshBlock permitting); F_cg(8 on C256) is the flux the C256
          flow carries at scales >= a C32 cell.
  --jobs  worker processes (login node: <= 16)
C256: point --rst/--bin at the C256 run, same --eost/--eosg as its metallicity, same
--nband, windows by restart index at the same rotations.

ENERGY (same as the run's etotgrav bookkeeping; src/pgen/deep_hot_jupiter_rt.cpp
TotPotAt, src/utils/atm_column.hpp GravPotAt): E = e + KE + rho*Phi,
  Phi = g ap (1 - ap/r) - 0.5 Omega^2 (r sin theta)^2  (point mass, rot_potential on)
The restart's conserved E already contains rho*Phi; e is the EINTRST1 slab; the script
checks E - e - KE - rho Phi(ours) = 0 (printed as phicheck).
Total specific enthalpy H = (E + p)/rho = h + k + Phi, h = (e + p)/rho, k = KE/rho.
Flux of the fluid through radial face f (between cells i-1 and i) at time t, column c:
  m_f = (m1_{i-1} + m1_i)/2 (m1 = rho v_r, the stored radial momentum),
  H_f = (H_{i-1} + H_i)/2  (same for h, k, Phi), F(t) = r_f^2 sum_c dOmega_c m_f H_f.
This is the RESOLVED advective flux; the Riemann (ledger FXE) flux also contains the
solver dissipation, so ledger - resolved = the numerical (dissipative) energy flux.

DECOMPOSITION (exact, per face, per energy component; no density weighting: m is the
mass flux, H the specific enthalpy, overbars are plain time means over the restarts of
the window at fixed column, [ ] the dOmega-weighted mean over a latitude band, * the
deviation from it, ' the deviation from the time mean at fixed column):
  The net-mass-flux part is removed FIRST, per restart: H~ = H - <H>(t, f) (<.> =
  dOmega-weighted mean over the face), NET = overline(Mdot(t) <H>(t)).  The C32 cell-
  centred m1 carries a large radially odd-even l=0 pattern (|Mdot| ~ 1e19 g/s, sign
  alternating face to face) whose NET = +-1e31 erg/s swamps everything; NET carries no
  heat (it depends on the zero of H) and is reported separately only.
  TOTAL = r^2 sum_c dOmega overline(m H~)    (the heat flux = total - NET)
  MMC   = r^2 sum_b Omega_b [mbar]_b [H~bar]_b                 (mean meridional)
  SE    = r^2 sum_b Omega_b [mbar* Hbar*]_b                     (stationary eddy)
  TR    = r^2 sum_c dOmega overline(m' H')  = TOTAL - MMC - SE  (transient, incl. the
          transients of the zonal means)
SCALES: F_cg(n) = the same TOTAL computed after replacing m and H by their dOmega-
weighted means over n x n cell blocks (within a MeshBlock, per restart); the flux carried
by scales below n cells is TOTAL - F_cg(n).
K_eff: F = -K 4 pi r^2 <rho><T> ds/dr (down-gradient entropy mixing; s per unit mass
from the EOS grid, shell means <.> = dOmega-weighted, time-averaged; ds/dr between the
two cell centroids).  K > 0 = down-gradient (what numerical mixing does), K < 0 =
counter-gradient.  oe = shell rms(v_r(i+1) - 2 v_r(i) + v_r(i-1))/4 / shell rms(v_r)
(1 = pure radial odd-even (2 dr) mode, ~0.61 = white noise, ~0 = smooth).  Grid scales: c_s dr, c_s dx_h, v_rms,r dr, |v|_rms dx_h with
dx_h = r (pi/2)/N (N = cells per panel edge), c_s = sqrt(Gamma1 p/rho).
"""
import argparse
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/ptmp2/jinma/deepconv_0925')
import dhjcs        # noqa: E402
import dhj_remap as R  # noqa: E402

RU = 8.314462618e7
BAR = 1e6
G = {}


def setup(a):
    h = R.read_rst(a.rst % a.first)
    par = h['par']
    nx1, r0, r1, c = R.grid_of(par)
    e = R.edges(nx1, r0, r1, c)
    rc = R.centroids(e)
    ng = h['ng']
    raw = dhjcs.bin_convert.read_binary(a.bin)
    geo = dhjcs.cs_geometry(raw)
    Gm = np.asarray(raw['mb_geometry'])
    nmb, nk, nj = h['nmb'], h['nx3'], h['nx2']
    dom = np.empty((nmb, nk, nj))
    cosc = np.empty((nmb, nk, nj))
    for m in range(nmb):
        x = np.tan(np.pi/4*(Gm[m, 2] + (Gm[m, 3] - Gm[m, 2])*(np.arange(nj) + 0.5)/nj))
        y = np.tan(np.pi/4*(Gm[m, 4] + (Gm[m, 5] - Gm[m, 4])*(np.arange(nk) + 0.5)/nk))
        X, Y = np.meshgrid(x, y)
        dxi = np.pi/4*(Gm[m, 3] - Gm[m, 2])/nj
        deta = np.pi/4*(Gm[m, 5] - Gm[m, 4])/nk
        dom[m] = dxi*deta*(1 + X*X)*(1 + Y*Y)/(1 + X*X + Y*Y)**1.5
        cosc[m] = -X*Y/np.sqrt((1 + X*X)*(1 + Y*Y))
    pb = par['problem']
    ap = float(pb['ap'])
    grav = float(pb['grav'])
    om = float(pb['omega'])
    rotpot = pb.get('rot_potential', 'false').lower() in ('true', '1')
    lat = geo['lat']
    phi = grav*ap*(1 - ap/rc)[None, None, None, :] \
        - 0.5*(om*om if rotpot else 0.0)*(rc[None, None, None, :]*np.cos(lat)[..., None])**2
    band = np.clip(((np.degrees(lat) + 90)//(180.0/a.nband)).astype(int), 0, a.nband - 1)
    npan = int(par['mesh']['nx2'])
    G.update(ng=ng, e=e, rc=rc, dom=dom, cosc=cosc, phi=phi, band=band, nx1=nx1,
             npan=npan, cg=[int(v) for v in a.cg.split(',')], nband=a.nband,
             eost=a.eost, eosg=a.eosg, rst=a.rst, prot=a.prot)
    return h['par']


def cgmean(q, w, n):
    """dOmega-weighted mean over n x n blocks of the (nmb, nk, nj, ...) arrays"""
    if n == 1:
        return q, w
    s = q.shape
    ww = w.reshape(s[0], s[1]//n, n, s[2]//n, n)
    qq = q.reshape(s[0], s[1]//n, n, s[2]//n, n, *s[3:])
    wx = ww.reshape(ww.shape + (1,)*(q.ndim - 3))
    W = ww.sum((2, 4))
    return (qq*wx).sum((2, 4))/W.reshape(W.shape + (1,)*(q.ndim - 3)), W


_eos = {}


def one(idx):
    h = R.read_rst(G['rst'] % idx)
    ng = G['ng']
    sl = (slice(None), slice(ng, -ng), slice(ng, -ng), slice(ng, -ng))
    u = h['u']
    rho = u[:, 0][sl[0], sl[1], sl[2], sl[3]]
    m1, m2, m3, E = [u[:, k][sl[0], sl[1], sl[2], sl[3]] for k in (1, 2, 3, 4)]
    p = h['p'][sl[0], sl[1], sl[2], sl[3]]
    ei = h['eint'][sl[0], sl[1], sl[2], sl[3]]
    c = G['cosc'][..., None]
    KE = R.kinetic(rho, m1, m2, m3, c)
    phi = G['phi']
    res = E - ei - KE - rho*phi
    phichk = np.abs(res).max()/np.abs(E).max()
    comp = {'h': (ei + p)/rho, 'k': KE/rho, 'phi': np.broadcast_to(phi, rho.shape)}
    f = lambda q: 0.5*(q[..., :-1] + q[..., 1:])   # noqa: E731  cell -> interior face
    mf = f(m1)                                      # (nmb, nk, nj, nx1-1)
    out = {'t': h['t'], 'phichk': phichk, 'm': mf}
    dom = G['dom']
    rf2 = G['e'][1:-1]**2
    Om = dom.sum()
    out['Mdot'] = rf2*(dom[..., None]*mf).sum((0, 1, 2))
    for k, q in comp.items():
        qf = f(q)
        gm = (dom[..., None]*qf).sum((0, 1, 2))/Om      # <H_k>(t, f)
        out['NET_' + k] = out['Mdot']*gm
        qf = qf - gm                                       # H~ = H - <H>(t, f)
        out['H_' + k] = qf
        out['mH_' + k] = mf*qf
    # coarse-grained totals (all components summed)
    Hf = out['H_h'] + out['H_k'] + out['H_phi']
    cg = []
    for n in G['cg']:
        mc, W = cgmean(mf, dom, n)
        Hc, _ = cgmean(Hf, dom, n)
        Hc = Hc - (W[..., None]*Hc).sum((0, 1, 2))/W.sum()
        cg.append(rf2*(W[..., None]*mc*Hc).sum((0, 1, 2)))
    out['cg'] = np.array(cg)
    # shell profiles for K_eff / grid scales (dOmega-weighted means per cell shell)
    if 'e' not in _eos:
        _eos['e'] = dhjcs.EOS(G['eost'])
        import dclib
        _eos['g'] = dclib.EOS(G['eosg'])
    T, _ = _eos['e'].invert(rho, ei)
    eg = _eos['g']
    s = eg.at('s', rho, T)*RU
    gad = eg.at('gad', rho, T)
    g1 = eg.at('chr', rho, T)/(1 - gad*eg.at('cht', rho, T))
    wn = (dom/dom.sum())[..., None]
    v2 = 2*KE/rho
    vr = m1/rho
    prof = {'rho': rho, 'T': T, 's': s, 'p': p, 'lnp': np.log(p),
            'cs2': g1*p/rho, 'vr2': vr*vr, 'v2': v2}
    out['prof'] = {k: (q*wn).sum((0, 1, 2)) for k, q in prof.items()}
    # radial 2-dr (odd-even) indicator: shell rms of the second difference of v_r / 4
    d2 = np.full(rho.shape[-1], np.nan)
    d2[1:-1] = ((vr[..., 2:] - 2*vr[..., 1:-1] + vr[..., :-2])**2*wn).sum((0, 1, 2))
    out['prof']['d2'] = d2/16.0
    out['prof']['mdotc'] = G['rc']**2*(dom[..., None]*m1).sum((0, 1, 2))
    return out


def rtime(fn):
    return R.read_rst(fn)['t']


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument('--rst', required=True)
    ap_.add_argument('--bin', required=True)
    ap_.add_argument('--eost', required=True)
    ap_.add_argument('--eosg', required=True)
    ap_.add_argument('--win', required=True)
    ap_.add_argument('--prot', type=float, default=1.101535e5)
    ap_.add_argument('--out', required=True)
    ap_.add_argument('--label', required=True)
    ap_.add_argument('--ledger', default=None)
    ap_.add_argument('--nband', type=int, default=30)
    ap_.add_argument('--cg', default='1,2,4,8')
    ap_.add_argument('--jobs', type=int, default=8)
    a = ap_.parse_args()
    wins = []
    for w in a.win.split(','):
        v = [int(x) for x in w.split(':')]
        wins.append(list(range(v[0], v[1] + 1, v[2] if len(v) > 2 else 1)))
    a.first = wins[0][0]
    setup(a)
    dom, band, nb = G['dom'], G['band'], G['nband']
    e = G['e']
    rf = e[1:-1]
    rf2 = rf**2
    # ledger face fluxes
    led = None
    if a.ledger:
        rows = [ln.split() for ln in open(a.ledger)
                if ln[:4].strip().isdigit() and len(ln.split()) == 13]
        if not rows or int(rows[0][0]) != 0:
            raise SystemExit('ledger has no complete per-shell table (needs rows from 0)')
        dfl = np.array([float(r[3]) for r in rows])
        led = -np.cumsum(dfl)            # FXE at face i+1 (FXE(0) = 0, closed wall)
    txt = []
    res = {}
    with Pool(a.jobs) as pool:
        for iw, idxs in enumerate(wins):
            outs = pool.map(one, idxs)
            nt = len(outs)
            rot = np.array([o['t'] for o in outs])/a.prot
            comps = ['h', 'k', 'phi']
            mbar = sum(o['m'] for o in outs)/nt
            TOT, MMC, NET, SE, TRs = {}, {}, {}, {}, {}
            ser = {}
            for k in comps:
                Hbar = sum(o['H_' + k] for o in outs)/nt
                mHbar = sum(o['mH_' + k] for o in outs)/nt
                TOT[k] = rf2*(dom[..., None]*mHbar).sum((0, 1, 2))
                ser[k] = np.array([rf2*(dom[..., None]*o['mH_' + k]).sum((0, 1, 2))
                                   for o in outs])
                mmc = np.zeros_like(rf)
                se = np.zeros_like(rf)
                for b in range(nb):
                    mk = band == b
                    if not mk.any():
                        continue
                    wb = dom[mk][:, None]
                    Ob = wb.sum()
                    mz = (wb*mbar[mk]).sum(0)/Ob
                    Hz = (wb*Hbar[mk]).sum(0)/Ob
                    mmc += rf2*Ob*mz*Hz
                    se += rf2*(wb*(mbar[mk] - mz)*(Hbar[mk] - Hz)).sum(0)
                MMC[k], SE[k] = mmc, se
                NET[k] = sum(o['NET_' + k] for o in outs)/nt
                TRs[k] = TOT[k] - mmc - se
            cg = sum(o['cg'] for o in outs)/nt
            prof = {q: sum(o['prof'][q] for o in outs)/nt for q in outs[0]['prof']}
            pf = np.exp(0.5*(prof['lnp'][:-1] + prof['lnp'][1:]))/BAR
            rc = G['rc']
            dsdr = np.diff(prof['s'])/np.diff(rc)
            rhoT = 0.5*(prof['rho'][:-1]*prof['T'][:-1] + prof['rho'][1:]*prof['T'][1:])
            D = 4*np.pi*rf2*rhoT*dsdr          # F = -K D
            tot = sum(TOT.values())
            eddy = sum(SE.values()) + sum(TRs.values())
            cs = np.sqrt(0.5*(prof['cs2'][:-1] + prof['cs2'][1:]))
            vr = np.sqrt(0.5*(prof['vr2'][:-1] + prof['vr2'][1:]))
            vv = np.sqrt(0.5*(prof['v2'][:-1] + prof['v2'][1:]))
            dr = np.diff(rc)
            dxh = rf*np.pi/2/G['npan']
            Ktot, Ked = -tot/D, -eddy/D
            oe_c = np.sqrt(prof['d2']/prof['vr2'])
            oe = 0.5*(oe_c[:-1] + oe_c[1:])      # 1 = pure 2-dr mode, 0.61 = white noise
            sert = sum(ser.values())
            sem = sert.std(0)/np.sqrt(nt)
            key = 'w%d' % iw
            res[key] = dict(rot=rot, pf=pf, rf=rf, TOT=TOT, MMC=MMC, NET=NET, SE=SE, TR=TRs,
                            cg=cg, cgn=np.array(G['cg']), series=ser, prof=prof, D=D,
                            Ktot=Ktot, Ked=Ked, cs=cs, oe=oe, vr=vr, vv=vv, dr=dr, dxh=dxh,
                            phichk=max(o['phichk'] for o in outs))
            txt.append('# window %d: %d restarts, rot %.2f..%.2f; phicheck max |E-e-KE-rho '
                       'Phi|/max|E| = %.1e' % (iw, nt, rot[0], rot[-1], res[key]['phichk']))
            txt.append('# erg/s, outward +; face f between cells f-1 and f; 4 pi factor: '
                       'fluxes are r^2 sum dOmega (full sphere)')
            txt.append('#  f  p[bar]  TOTAL(+-sem)  | MMC(NETnotincl) SE TR | h k phi (TOTAL) |'
                       ' %s | ledgerFXE ledger-TOTAL oe | K_tot K_eddy [cm2/s] | K/(cs dr) '
                       'K/(cs dxh) K/(vr dr) K/(|v| dxh)'
                       % ' '.join('cg%d' % n for n in G['cg'][1:]))
            for i in range(len(rf)):
                sm = lambda d_: sum(d_[k][i] for k in comps)   # noqa: E731
                lg = led[i] if led is not None and i < len(led) else np.nan
                cgs = ' '.join('%+.2e' % cg[j, i] for j in range(1, len(G['cg'])))
                txt.append('%3d %9.3g %+.3e(%.1e) | %+.3e(%+.2e) %+.3e %+.3e | %+.3e %+.3e '
                           '%+.3e | %s | %+.3e %+.2e %.2f | %+.2e %+.2e | %+.2e %+.2e %+.2e %+.2e' % (
                               i + 1, pf[i], tot[i], sem[i], sm(MMC), sm(NET), sm(SE),
                               sm(TRs), TOT['h'][i], TOT['k'][i], TOT['phi'][i], cgs, lg, lg - tot[i], oe[i],
                               Ktot[i], Ked[i], Ktot[i]/(cs[i]*dr[i]),
                               Ktot[i]/(cs[i]*dxh[i]), Ktot[i]/(vr[i]*dr[i]),
                               Ktot[i]/(vv[i]*dxh[i])))
            print('\n'.join(txt[-(len(rf) + 3):]), flush=True)
    np.savez(os.path.join(a.out, 'deepmix_%s.npz' % a.label),
             res=np.array(res, dtype=object), ledger=led)
    open(os.path.join(a.out, 'deepmix_%s.txt' % a.label), 'w').write('\n'.join(txt) + '\n')


if __name__ == '__main__':
    main()
