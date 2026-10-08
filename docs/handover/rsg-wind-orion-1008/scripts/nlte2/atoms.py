"""Statistical-equilibrium model atoms (Kurucz gfemq E1+M1+E2 line lists) in the diluted
stellar field, Sobolev escape, collisions with H/H2 (Drawin x S_H, Barklem 2018 for Fe I,
forbidden-manifold quench q_forb, fine structure q_fs) and electrons (van Regemorter,
Omega=1 forbidden).  Returns the net collisional heating per atom (erg/s).

Sign: heating = sum_pairs (n_u C_ul - n_l C_lu)(E_u - E_l)   (> 0 heats the gas).
"""
import re
import numpy as np

D = '/orion/ptmp/jinma/rsg_wind_1008/nlte2/data/'
H = 6.62607e-27
C = 2.99792458e10
KB = 1.380649e-16
EV = 1.602177e-12
HCK = 1.438777       # cm K
MH = 1.6726e-24
ME = 9.109e-28
A0 = 0.529177e-8
EH = 13.6057 * EV
PIE2MC = 0.026540    # pi e^2 / m c  (cm^2 Hz)
LC = 'SPDFGHIKLMNOQ'

# element data: mass (amu), IP (eV), ground-state photoionisation threshold cross section
# (cm^2; rounded Verner+1996 / TOPbase order-of-magnitude values, uncertain x3-10)
ELEM = {
    'Fe': (55.85, 7.902, 5e-19), 'Ca': (40.08, 6.113, 4e-19), 'Ti': (47.87, 6.828, 1e-18),
    'Cr': (52.00, 6.767, 5e-18), 'Mn': (54.94, 7.434, 1e-18), 'Ni': (58.69, 7.640, 1e-18),
    'Mg': (24.31, 7.646, 1.2e-18), 'Na': (22.99, 5.139, 1.1e-19), 'K': (39.10, 4.341, 1e-20),
    'Al': (26.98, 5.986, 6.5e-17), 'Si': (28.09, 8.152, 3.7e-17),
    'O': (16.00, 13.618, 0.), 'C': (12.01, 11.260, 0.),
}
CODE = {'Fe': 26, 'Ca': 20, 'Ti': 22, 'Cr': 24, 'Mn': 25, 'Ni': 28, 'Mg': 12, 'Na': 11,
        'K': 19, 'Al': 13, 'Si': 14, 'O': 8, 'C': 6}


def term_of(label):
    """'4s2 a5D' -> ('a', 5, 'D');  '(2F)4p 3F' -> ('', 3, 'F')."""
    tok = label.split()[-1] if label.split() else ''
    m = re.match(r'^([a-z]?)(\d)([A-Z])', tok)
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


class Atom:
    def __init__(self, el, ion=0, emax_ev=6.0, kind='pos'):
        self.el, self.ion = el, ion
        z = CODE[el]
        fn = D + f'gfemq{z:02d}{ion:02d}.pos'
        self.mass = ELEM[el][0]
        emax = emax_ev * EV / (H * C)
        lev = {}
        raw = []
        for s in open(fn):
            try:
                e1, j1, e2, j2 = float(s[24:36]), float(s[36:41]), float(s[52:64]), float(s[64:69])
                lgf = float(s[11:18])
            except ValueError:
                continue
            e1, e2 = abs(e1), abs(e2)
            if max(e1, e2) > emax:
                continue
            k1, k2 = (round(e1, 1), j1), (round(e2, 1), j2)
            for k, e, j, lab in ((k1, e1, j1, s[42:52]), (k2, e2, j2, s[70:80])):
                if k not in lev:
                    lev[k] = (e, j, lab.strip())
            raw.append((k1, k2, lgf, s[101]))
        keys = sorted(lev, key=lambda k: lev[k][0])
        idx = {k: i for i, k in enumerate(keys)}
        self.E = np.array([lev[k][0] for k in keys])            # cm^-1
        self.g = np.array([2*lev[k][1] + 1 for k in keys])
        self.lab = [lev[k][2] for k in keys]
        self.N = len(keys)
        # merge lines per level pair
        L = {}
        for k1, k2, lgf, typ in raw:
            i, j = idx[k1], idx[k2]
            lo, up = (i, j) if self.E[i] < self.E[j] else (j, i)
            if lo == up or self.E[up] - self.E[lo] < 1e-3:
                continue
            gf = 10**lgf
            sig = self.E[up] - self.E[lo]
            A = 0.6670e16 * gf / (self.g[up] * (1e8/sig)**2)
            d = L.setdefault((lo, up), [0., 0., False])
            d[0] += A
            if typ == ' ':
                d[1] += gf / self.g[lo]
                d[2] = True
        pr = np.array(list(L.keys()))
        self.il, self.iu = pr[:, 0], pr[:, 1]
        self.A = np.array([v[0] for v in L.values()])
        self.f = np.array([v[1] for v in L.values()])
        self.e1 = np.array([v[2] for v in L.values()])
        self.sig = self.E[self.iu] - self.E[self.il]
        self.lam = 1/self.sig                                   # cm
        # parity by 2-colouring the E1 graph from the ground level
        par = -np.ones(self.N, int)
        par[0] = 0
        adj = [[] for _ in range(self.N)]
        for a, b in zip(self.il[self.e1], self.iu[self.e1]):
            adj[a].append(b)
            adj[b].append(a)
        stack = [0]
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if par[b] < 0:
                    par[b] = 1 - par[a]
                    stack.append(b)
        self.par = par
        # keep only levels connected to the ground by E1 (others: isolated, drop)
        keep = par >= 0
        if not keep.all():
            self._subset(keep)
        eodd = self.E[self.par != self.par[0]].min() if (self.par != self.par[0]).any() else 0.
        adn = np.zeros(self.N)
        np.add.at(adn, self.iu[self.e1], self.A[self.e1])
        # metastable = ground parity and total E1 downward A < 1e3 /s (low even manifold)
        self.meta = (self.par == self.par[0]) & (adn < 1e3)
        self.eodd = eodd
        self.term = [term_of(x) for x in self.lab]
        self._pairs()

    def _subset(self, keep):
        new = -np.ones(self.N, int)
        new[keep] = np.arange(keep.sum())
        m = keep[self.il] & keep[self.iu]
        for a in ('il', 'iu'):
            setattr(self, a, new[getattr(self, a)[m]])
        for a in ('A', 'f', 'e1', 'sig', 'lam'):
            setattr(self, a, getattr(self, a)[m])
        self.E, self.g, self.par = self.E[keep], self.g[keep], self.par[keep]
        self.lab = [x for x, k in zip(self.lab, keep) if k]
        self.N = keep.sum()

    def _pairs(self):
        """collision pair list: E1 lines (Drawin / van Regemorter) + forbidden pairs inside
        the low even manifold (q_forb, fine structure q_fs if same term)."""
        pl, pu, fe1, kind = [], [], [], []
        for a, b, f, e in zip(self.il, self.iu, self.f, self.e1):
            if e and f > 0:
                pl.append(a), pu.append(b), fe1.append(f), kind.append(0)
        mi = np.where(self.meta)[0]
        for x in range(len(mi)):
            for y in range(x+1, len(mi)):
                a, b = mi[x], mi[y]
                ta, tb = self.term[a], self.term[b]
                same = (ta is not None and ta == tb and
                        abs(self.E[a] - self.E[b]) < 2500.)
                pl.append(a), pu.append(b), fe1.append(0.), kind.append(2 if same else 1)
        self.pl, self.pu = np.array(pl), np.array(pu)
        self.pf, self.pk = np.array(fe1), np.array(kind)
        self.pdE = (self.E[self.pu] - self.E[self.pl]) * H * C       # erg
        self.bark = None

    def add_barklem(self):
        """Fe I: Barklem (2018) term-resolved H rates; level i of term I gets R_IJ g_j/g_J.
        Replaces Drawin for the pairs where it is non-zero."""
        b = D + 'barklem18/'
        S = [s.split() for s in open(b + 'states.dat')][:166]
        Eb = np.array([float(s[-1]) for s in S])
        lb = [s[1] for s in S]
        tb = []
        for x in lb:
            m = re.match(r'^(?:\w+\.)?([a-z]?)(\d)([A-Z])', x)
            tb.append((int(m.group(2)), m.group(3), x.rstrip('0123456789').endswith('o')))
        ip = (self.par != self.par[0])
        mapl = -np.ones(self.N, int)
        for i in range(self.N):
            t = self.term[i]
            if t is None:
                continue
            best, bd = -1, 1500.
            for k in range(166):
                if tb[k][0] == t[1] and tb[k][1] == t[2] and tb[k][2] == bool(ip[i]):
                    d = abs(Eb[k] - self.E[i])
                    if d < bd:
                        best, bd = k, d
            mapl[i] = best
        gT = np.zeros(166)
        np.add.at(gT, mapl[mapl >= 0], self.g[mapl >= 0])
        Ts = np.array([1000, 2000, 3000, 4000, 5000.])
        Rt = np.array([np.loadtxt(b + f'{int(t)}_K.rates')[:166, :166] for t in Ts])
        self.bmap, self.bgT, self.bTs, self.bR = mapl, gT, Ts, Rt
        self.nbark = (mapl >= 0).sum()
        # add Barklem-only term pairs (non-E1, non-manifold) to the pair list
        have = set(zip(self.pl.tolist(), self.pu.tolist()))
        nz = (Rt[2] > 0)
        add_l, add_u = [], []
        lv = [np.where(mapl == k)[0] for k in range(166)]
        for I, J in zip(*np.nonzero(nz)):
            for i in lv[J]:            # initial term J -> final term I
                for j in lv[I]:
                    a, c = (i, j) if self.E[i] < self.E[j] else (j, i)
                    if a != c and (a, c) not in have:
                        have.add((a, c))
                        add_l.append(a), add_u.append(c)
        if add_l:
            self.pl = np.concatenate([self.pl, add_l])
            self.pu = np.concatenate([self.pu, add_u])
            self.pf = np.concatenate([self.pf, np.zeros(len(add_l))])
            self.pk = np.concatenate([self.pk, np.full(len(add_l), 3)])
            self.pdE = (self.E[self.pu] - self.E[self.pl]) * H * C
        self.bark = True

    def barklem_down(self, T):
        """downward H rate coefficient per pair (cm^3/s), 0 where none."""
        lt = np.log(np.clip(T, 1000, 5000))
        lTs = np.log(self.bTs)
        k = np.clip(np.searchsorted(lTs, lt) - 1, 0, 3)
        w = (lt - lTs[k]) / (lTs[k+1] - lTs[k])
        I = self.bmap[self.pu]
        Jl = self.bmap[self.pl]
        ok = (I >= 0) & (Jl >= 0) & (I != Jl)
        q = np.zeros(len(self.pl))
        r0 = self.bR[k][Jl[ok], I[ok]]
        r1 = self.bR[k+1][Jl[ok], I[ok]]
        with np.errstate(divide='ignore'):
            lr = (1-w)*np.log(np.maximum(r0, 1e-300)) + w*np.log(np.maximum(r1, 1e-300))
        q[ok] = np.where((r0 > 0) & (r1 > 0), np.exp(lr), 0.) * self.g[self.pl[ok]] / \
            np.maximum(self.bgT[Jl[ok]], 1)
        return q

    # ---------------------------------------------------------------- rates
    def coll(self, T, nH, nH2, ne, SH, qforb, qfs, fH2=1.0):
        """downward collision rates per pair (1/s) and upward via detailed balance."""
        x = self.pdE / (KB*T)
        mu = self.mass*MH*MH/(self.mass*MH + MH)
        # Drawin (Lambert 1993 form), upward rate coefficient
        dE = np.maximum(self.pdE, 1e-30)
        psi = np.exp(-x) / (1 + 2/x)
        qup_d = (16*np.pi*A0**2*np.sqrt(2*KB*T/(np.pi*mu))*(EH/dE)**2*self.pf
                 * self.mass*ME/(ME + MH) * psi)
        gl, gu = self.g[self.pl], self.g[self.pu]
        qdn_d = qup_d * gl/gu * np.exp(x)
        qdn_H = SH * qdn_d
        if self.bark:
            qb = self.barklem_down(T)
            qdn_H = np.where(qb > 0, qb, qdn_H)
        qdn_H = np.where(self.pk == 1, qforb, qdn_H)
        qdn_H = np.where(self.pk == 2, qfs, qdn_H)
        # electrons: van Regemorter for E1, Omega = 1 for forbidden
        from scipy.special import exp1
        gbar = 0.276*np.exp(np.minimum(x, 700))*exp1(np.minimum(x, 700))
        if self.ion > 0:
            gbar = np.maximum(gbar, 0.2)
        Om = 14.5*gl*self.pf*(EH/dE)*gbar
        Om = np.where(self.pk > 0, 1.0, Om)
        qdn_e = 8.629e-6*Om/(gu*np.sqrt(T))
        Cd = (nH + fH2*nH2)*qdn_H + ne*qdn_e
        Cu = Cd * gu/gl * np.exp(-x)
        return Cd, Cu

    def solve(self, T, nH, nH2, ne, nX, W, Teff, dvdr, SH=0.1, qforb=1e-11, qfs=1e-9,
              uvfac=1.0, fH2=1.0, niter=12, lte=False):
        """populations (normalised), heating per atom, absorbed power per atom, betas."""
        nu = C*self.sig
        x_s = HCK*self.sig/Teff
        nbar = W/np.expm1(x_s) * np.where(self.lam < 3e-5, uvfac, 1.0)
        Cd, Cu = self.coll(T, nH, nH2, ne, SH, qforb, qfs, fH2)
        gl, gu = self.g[self.il], self.g[self.iu]
        boltz = self.g*np.exp(-HCK*self.E/T)
        boltz /= boltz.sum()
        beta = np.ones(len(self.A))
        n = boltz.copy()
        N = self.N
        for it in range(niter if dvdr > 0 else 1):
            if lte:
                n = boltz
            else:
                M = np.zeros((N, N))
                Rd = self.A*(1 + nbar)*beta
                Ru = self.A*nbar*beta*gu/gl
                np.add.at(M, (self.iu, self.il), Rd)        # M[i,j] = rate i->j
                np.add.at(M, (self.il, self.iu), Ru)
                np.add.at(M, (self.pu, self.pl), Cd)
                np.add.at(M, (self.pl, self.pu), Cu)
                Q = M.T - np.diag(M.sum(1))
                Q[0, :] = 1.
                rhs = np.zeros(N)
                rhs[0] = 1.
                n = np.linalg.solve(Q, rhs)
                n = np.maximum(n, 0.)
            if dvdr <= 0:
                break
            nl, nu_ = n[self.il]*nX, n[self.iu]*nX
            tau = self.lam**3*self.A*gu/(8*np.pi*gl)*np.maximum(nl - nu_*gl/gu, 0.)/dvdr
            bnew = np.where(tau > 1e-6, -np.expm1(-tau)/np.maximum(tau, 1e-300), 1.)
            if np.allclose(bnew, beta, rtol=1e-3, atol=1e-12) and it > 0:
                beta = bnew
                break
            beta = np.sqrt(beta*bnew) if it > 2 else bnew
        heat = ((n[self.pu]*Cd - n[self.pl]*Cu)*self.pdE).sum()
        hnu = H*nu
        absorb = (self.A*beta*hnu*(n[self.il]*gu/gl*nbar - n[self.iu]*nbar)).sum()
        emit = (self.A*beta*hnu*n[self.iu]).sum()
        return dict(n=n, heat=heat, absorb=absorb, emit=emit, beta=beta, nbar=nbar)

    def heat_lte(self, T, W, Teff, beta=None, uvfac=1.0):
        """LTE (Boltzmann at T over the model levels) line heating per atom, thin or x beta."""
        nu = C*self.sig
        x_s = HCK*self.sig/Teff
        x_t = HCK*self.sig/T
        nbar = W/np.expm1(x_s) * np.where(self.lam < 3e-5, uvfac, 1.0)
        b = self.g*np.exp(-HCK*self.E/T)
        b /= b.sum()
        gl, gu = self.g[self.il], self.g[self.iu]
        bb = 1. if beta is None else beta
        return (self.A*bb*H*nu*b[self.il]*gu/gl*(nbar - np.exp(-x_t)*(1 + nbar))).sum()

    def gamma_pi(self, n, W, Teff, uvfac=1.0, sig_exc=3e-18):
        """photoionisation rate per atom (1/s) in W B(Teff), sigma ~ sig_th (nu_th/nu)^2."""
        IP = ELEM[self.el][1]*EV/(H*C)
        sth = np.where(np.arange(self.N) == 0, ELEM[self.el][2], sig_exc)
        lam = np.geomspace(5e-6, 1.2e-4, 600)       # cm, 50-1200 nm
        sig = 1/lam
        nphot = 8*np.pi*C/lam**4/np.expm1(np.minimum(HCK/(lam*Teff), 700)) * W
        nphot *= np.where(lam < 3e-5, uvfac, 1.0)    # photons /cm^2/s per cm (4 pi J/h nu)
        g = 0.
        for i in np.where(n > 1e-12)[0]:
            th = IP - self.E[i]
            m = sig > th
            if m.sum() < 2:
                continue
            g += n[i]*sth[i]*np.trapezoid(nphot*(th/sig)**2*m, lam)
        return abs(g)
