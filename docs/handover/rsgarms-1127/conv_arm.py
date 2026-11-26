# conv_arm.py: RSGARMS copy of radtests/oldnew/conv_on.py; only change: the vet dump test matches any dump number (restarted NEW dirs have no 00000 dump)
"""conv.py KEY RUNDIR OUT.npz : shell profiles of every dump of a long star3d run (sp wedge or cubed sphere): mean state, N^2, velocities,
radiative / convective / kinetic luminosities.  EOS: sp gamma law 5/3 (code T = P/rho); cs the dumped table r_eos/rsg_eos.txt."""
import os, sys, glob, json, re
import numpy as np
key, D, out = sys.argv[1:4]; MV = 'vet' if glob.glob(D + '/bin/*.vet.*.bin') else 'm1'; HERE = os.path.dirname(os.path.abspath(__file__)); DAY = 86400.0; RSUN = 6.957e10
from concurrent.futures import ProcessPoolExecutor
CSG = key in ('cs', 'bsgcs'); TAB = key == 'cs'
C = dict(sp=dict(tu=7.4568864364e-9, gm=2.62502e27), bsgcs=dict(tu=7.4568864364e-9, gm=2.62502e27), cs=dict(tu=1.202724e-8, gm=1.4842265008e27))[key]
SH = dict(sp=(60, 78, 90, 96, 108, 126, 150), bsgcs=(60, 78, 90, 96, 108, 126, 150), cs=(110, 220, 264, 286, 308, 316, 330))[key]
ARK = 7.565723e-15; CL = 2.99792458e10
def read(path):
    with open(path, 'rb') as fp:
        size = os.fstat(fp.fileno()).st_size; fp.readline(); npre = int(fp.readline().split(b'=')[-1]); par = {}
        for _ in range(npre - 1):
            k, v = [x.strip() for x in fp.readline().decode().split('=')]; par[k] = v
        fp.readline(); names = fp.readline().decode().split()[1:]; fp.read(int(fp.readline().split(b'=')[-1])); bl = []
        while fp.tell() < size:
            ind = np.frombuffer(fp.read(24), np.int32); lx = np.frombuffer(fp.read(16), np.int32); fp.read(48)
            n1, n2, n3 = ind[1]-ind[0]+1, ind[3]-ind[2]+1, ind[5]-ind[4]+1; x = []
            for n in (n1, n2, n3): x.append((np.fromfile(fp, np.float64, n+1), np.fromfile(fp, np.float64, n)))
            a = np.fromfile(fp, np.float32 if int(par['size of variable']) == 4 else np.float64, n1*n2*n3*len(names)).astype(np.float64).reshape(len(names), n3, n2, n1)
            bl.append((lx[:3].copy(), x, dict(zip(names, a))))
    return float(par['time']), int(par['cycle']), bl
def cols(bl, names):
    """(6 or 1, N3, N2, n1) arrays + faces + solid-angle weights of the columns"""
    npan = 6 if CSG else 1; per = len(bl)//npan; n3, n2, n1 = bl[0][2][names[0]].shape
    L2 = max(b[0][1] for b in bl) + 1; L3 = max(b[0][2] for b in bl) + 1; o = {v: np.zeros((npan, L3*n3, L2*n2, n1)) for v in names}
    f2, f3 = np.zeros(L2*n2+1), np.zeros(L3*n3+1)
    for m, b in enumerate(bl):
        j, k = b[0][1]*n2, b[0][2]*n3; f2[j:j+n2+1] = b[1][1][0]; f3[k:k+n3+1] = b[1][2][0]
        for v in names: o[v][m//per, k:k+n3, j:j+n2] = b[2][v]
    if not CSG: W = np.outer(np.diff(f3), -np.diff(np.cos(f2)))[None]
    else:
        s = 1.0 if abs(f2[-1] - 1) < 0.3 else np.pi/4; s2, s3 = f2/ (4/np.pi) if abs(f2[-1]-1) < 1e-6 else f2, f3/(4/np.pi) if abs(f3[-1]-1) < 1e-6 else f3
        x, y = np.tan(0.5*(s2[1:]+s2[:-1])), np.tan(0.5*(s3[1:]+s3[:-1])); d = np.sqrt(1 + x[None]**2 + y[:, None]**2)
        W = (1+x[None]**2)*(1+y[:, None]**2)/d**3*np.diff(s2)[None]*np.diff(s3)[:, None]; W = np.broadcast_to(W, (6,)+W.shape)
    return o, bl[0][1][0][0], bl[0][1][0][1], W, (f2, f3)
class Eos:
    """cs: the run's own table; sp: gamma law.  P(rho,T), e(rho,T) [per gram] and T(rho,e)"""
    def __init__(s):
        if TAB:
            f = os.path.join(HERE, 'fpcpu/r_eos/rsg_eos.txt'); f = f if os.path.exists(f) else os.path.join(HERE, 'fpcpu/rsg_eos.txt')
            h = [l for l in open(f) if l.startswith('#')][2].split()[1:]; s.nx, s.ny = int(h[0]), int(h[1]); s.x0, s.dx, s.y0, s.dy = map(float, h[2:6])
            a = np.loadtxt(f); s.le, s.lp = a[:, 0].reshape(s.ny, s.nx), a[:, 1].reshape(s.ny, s.nx)
    def ip(s, tab, lr, lt):
        x = (lr - s.x0)/s.dx; y = (lt - s.y0)/s.dy; i = np.clip(x.astype(int), 0, s.nx-2); j = np.clip(y.astype(int), 0, s.ny-2); u, w = x - i, y - j
        return (1-w)*((1-u)*tab[j, i] + u*tab[j, i+1]) + w*((1-u)*tab[j+1, i] + u*tab[j+1, i+1])
    def T_of(s, rho, e):
        if not TAB: return 2/3*e/rho*C['tu']
        lr, le = np.log10(rho), np.log10(e/rho); lo, hi = np.full(rho.shape, s.y0), np.full(rho.shape, s.y0 + s.dy*(s.ny-1))
        for _ in range(45):
            m = 0.5*(lo + hi); up = s.ip(s.le, lr, m) < le; lo = np.where(up, m, lo); hi = np.where(up, hi, m)
        return 10**(0.5*(lo + hi))
    def P(s, rho, T): return rho*T/C['tu'] if not TAB else rho*10**s.ip(s.lp, np.log10(rho), np.log10(T))
    def e(s, rho, T): return 1.5*T/C['tu'] if not TAB else 10**s.ip(s.le, np.log10(rho), np.log10(T))
    def gad(s, rho, T, Er):
        """Gamma1, nabla_ad of gas + radiation in thermal equilibrium, and cp"""
        a3 = ARK; d = 1e-3
        Pt = lambda r, t: s.P(r, t) + ARK*t**4/3; et = lambda r, t: s.e(r, t) + ARK*t**4/r
        Pr, Ptl = Pt(rho*(1+d), T), Pt(rho, T*(1+d)); P0 = Pt(rho, T); chr_ = (np.log(Pr) - np.log(P0))/np.log(1+d); cht = (np.log(Ptl) - np.log(P0))/np.log(1+d)
        er = (et(rho*(1+d), T) - et(rho*(1-d), T))/(2*d)/rho; cv = (et(rho, T*(1+d)) - et(rho, T*(1-d)))/(2*d*T)
        dTdr = (P0/rho**2 - er)/cv*rho/T         # dlnT/dlnrho at constant s
        G1 = chr_ + cht*dTdr; nad = dTdr/G1; cp = cv + P0/(rho*T)*cht**2/chr_
        return G1, nad, cp, cht/chr_, P0
def prof(nm):
    t, cyc, bh = read(os.path.join(D, 'bin/%s.hydro_w.%s.bin' % (PRE, nm))); _, _, bv = read(os.path.join(D, 'bin/%s.%s.%s.bin' % (PRE, MV, nm)))
    H, xf, xc, W, ang = cols(bh, ['dens', 'velx', 'vely', 'velz', 'eint']); nE, nP = ('vet_e', 'vet_phi1') if MV == 'vet' else ('m1_e', 'm1_f1'); V, _, _, _, _ = cols(bv, [nE, nP]); V = dict(vet_e=V[nE], vet_phi1=V[nP])
    w = (W/W.sum())[..., None]; av = lambda q: (w*q).sum((0, 1, 2)); rho, v1, v2, v3, e = H['dens'], H['velx'], H['vely'], H['velz'], H['eint']; E, phi = V['vet_e'], V['vet_phi1']
    T = EOS.T_of(rho, e); P = EOS.P(rho, T)
    if CSG:
        s2, s3 = [f*np.pi/4 if abs(f[-1]-1) < 1e-6 else f for f in ang]; x, y = np.tan(0.5*(s2[1:]+s2[:-1])), np.tan(0.5*(s3[1:]+s3[:-1]))
        c = (-x[None]*y[:, None]/np.sqrt((1+x[None]**2)*(1+y[:, None]**2)))[None, :, :, None]; vl2 = v2**2 + v3**2 + 2*c*v2*v3
    else: vl2 = v2**2 + v3**2
    m, mv = av(rho), av(rho*v1); vm = mv/m; d = {}
    d['rho'], d['T'], d['Tmean_m'], d['P'], d['E'] = m, av(T), av(rho*T)/m, av(P), av(E); d['vmean'] = vm
    d['vr_rms'] = np.sqrt(av(rho*(v1 - vm)**2)/m); d['vl_rms'] = np.sqrt(av(rho*vl2)/m); d['dT_rms'] = np.sqrt(av((T/d['T'] - 1)**2)); d['drho_rms'] = np.sqrt(av((rho/m - 1)**2))
    h = e + P; hr = 4/3*E                                    # enthalpy densities (gas e + P; radiation E + P)
    d['Fg'] = av(h*v1) - av(h)*vm; d['Frad_adv'] = av(hr*v1) - av(hr)*vm        # convective: fluctuating (relative to the mean mass flow)
    d['Fk'] = av(0.5*rho*((v1)**2 + vl2)*v1); d['Fph'] = av(phi); d['Fdiff'] = av(phi - hr*v1); d['Fmean'] = (av(h) + av(hr))*vm
    d['ke_r'], d['ke_l'] = av(0.5*rho*(v1 - vm)**2), av(0.5*rho*vl2); d['vmax'] = np.abs(v1).max((0, 1, 2)); d['up'] = (w*(v1 > vm)).sum((0, 1, 2))
    vp = (v1 - vm).astype(np.float32); tp = (T/d['T'] - 1).astype(np.float32)
    if not CSG:                  # meridional (mid phi) and equatorial (mid theta) planes, shells
        n = v1.shape[1]//2; d['sv'], d['st'] = vp[0, n], tp[0, n]; d['ev'], d['et'] = vp[0, :, n], tp[0, :, n]
        d['shv'] = np.array([vp[0, :, :, i] for i in SH]); d['sht'] = np.array([tp[0, :, :, i] for i in SH])
    else:                            # equatorial ring: panels 0 1 2 4, the two rows about eta = 0
        n = v1.shape[1]//2; d['sv'] = np.concatenate([vp[p, n-1:n+1].mean(0) for p in (0, 1, 2, 4)]); d['st'] = np.concatenate([tp[p, n-1:n+1].mean(0) for p in (0, 1, 2, 4)])
        d['shv'] = np.array([vp[..., i] for i in SH]); d['sht'] = np.array([tp[..., i] for i in SH])
    d['t'], d['cyc'] = t, cyc; return d
if __name__ == '__main__':
    PRE = {'sp': 'bsg3d', 'cs': 'rsg', 'bsgcs': 'bsgcs'}[key]; EOS = Eos()
    nms = sorted(re.search(r'\.(\d{5})\.bin', f).group(1) for f in glob.glob(os.path.join(D, 'bin/%s.hydro_w.*.bin' % PRE)) if os.path.exists(f.replace('hydro_w', MV)))
    _, _, b0 = read(os.path.join(D, 'bin/%s.hydro_w.%s.bin' % (PRE, nms[0]))); xf, xc = b0[0][1][0]
    with ProcessPoolExecutor(int(os.environ.get('NP', '8'))) as ex: L = list(ex.map(prof, nms))
    S = {k: np.array([x[k] for x in L]) for k in L[0]}; S['xf'], S['r'] = xf, xc; g = C['gm']/xc**2; np.seterr(all='ignore')
    G = [EOS.gad(S['rho'][i], S['Tmean_m'][i], S['E'][i]) for i in range(len(L))]; S['G1'], S['nad'], S['cp'], S['Q'], _ = [np.array([q[j] for q in G]) for j in range(5)]
    Pt = S['P'] + S['E']/3; dl = lambda q: np.gradient(np.log(q), xc, axis=1); S['Ptot'] = Pt; S['nab'] = dl(S['T'])/dl(Pt); S['N2'] = g*(dl(Pt)/S['G1'] - dl(S['rho'])); S['Hp'] = -1/dl(Pt)
    S['krho'] = -CL/3*np.gradient(S['E'], xc, axis=1)/np.where(S['Fdiff'] > 0, S['Fdiff'], np.nan)     # effective kappa*rho from F = -c dE/dr/(3 kappa rho)
    S['W'] = cols(b0, ['dens'])[3]; S['th'] = np.concatenate([[0.0]]); S['SH'] = np.array(SH); S['ang2'], S['ang3'] = cols(b0, ['dens'])[4]; np.savez_compressed(out, **S); print(len(L), 'dumps', S['t'][[0, -1]]/DAY)
