"""ck11.py (synth2_emis): the GCM's own 11-band Exo-FMS opacity, copied from
/viper/ptmp2/jinma/wasp121_0925/rr_lib_met.py (lines 1-210) with explicit table paths (hiT2).

Original docstring: Real-ray check of the correlated-k stellar beam: shared library.

Python transcription of the model's own shortwave opacity (src/utils/correlated_k.hpp:
read_ck_table, read_ck_continuum, ck_tp_index, ck_kappa, ck_continuum, the Rosseland
table of ck_build_rosseland_table) and of the two beam rules in
src/utils/two_stream_rt.hpp (plane-parallel facsw path and the ck_beam_sph block).
Read-only on every input.
"""
import numpy as np

REPO = '/viper/u2/jinma/ATHENAK/athenak'
DATA = REPO + '/data/exo_fms_ck'
import os as _os
MET = _os.environ.get('RCE_MET', '1x')    # wasp121_0925: premixed metallicity tag
CKDIR = _os.environ.get('RCE_CKDIR', DATA)   # where ck/ and CE_tables/ for MET live
NB, NG = 11, 8
KB = 1.380649e-16
SIG = 5.6704e-5


def planck_below(lamT):
    lamT = np.atleast_1d(np.asarray(lamT, float))
    out = np.zeros_like(lamT)
    for q, x in enumerate(lamT):
        if x <= 0:
            continue
        xi = 1.4387769e4/x
        if xi > 700:
            continue
        s = 0.0
        for n in range(1, 501):
            nx = n*xi
            if nx > 700:
                break
            rn = 1.0/n
            s += np.exp(-nx)*rn*(xi**3 + 3*xi*xi*rn + 6*xi*rn*rn + 6*rn**3)
        out[q] = 15.0/np.pi**4*s
    return out


def tp_index(lg, x):
    """ck_tp_index, vectorised: clamped bracket."""
    x = np.asarray(x, float)
    n = len(lg)
    i = np.clip(np.searchsorted(lg, x, side='right') - 1, 0, n-2)
    f = (x - lg[i])/(lg[i+1] - lg[i])
    lo = ~(x > lg[0])
    hi = x >= lg[-1]
    i = np.where(lo, 0, np.where(hi, n-2, i))
    f = np.where(lo, 0.0, np.where(hi, 1.0, f))
    return i, f


class CK:
    def __init__(self, star_teff=6000.0, ktab=None, cetab=None, ddir=None):
        global DATA
        if ddir is not None:
            DATA = ddir
        with open(ktab) as fh:
            fh.readline()
            v = np.array(fh.read().split(), float)
        nT, nP, nb, ng = v[:4].astype(int)
        assert nb == NB and ng == NG
        p = 4
        self.lT = np.log10(v[p:p+nT]); p += nT
        self.lP = np.log10(v[p:p+nP]); p += nP
        self.wl = v[p:p+NB+1]; p += NB+1
        p += NB+1
        p += NG
        self.gw = v[p:p+NG]; p += NG
        k = v[p:p+nT*nP*NB*NG].reshape(nT, nP, NB, NG)
        assert k.size == nT*nP*NB*NG
        self.lk = np.log10(np.maximum(k, 1e-99)).transpose(2, 3, 0, 1).copy()  # b,g,T,P
        # FastChem
        with open(cetab) as fh:
            hdr = fh.readline().split()
            fh.readline()
            w = np.array(fh.read().split(), float)
        cT, cP = int(hdr[0]), int(hdr[1])
        self.ceT = np.log10(w[:cT]); self.ceP = np.log10(w[cT:cT+cP])
        self.ce = w[cT+cP:cT+cP+cT*cP*6].reshape(cT, cP, 6)
        # CIA
        self.cia = []
        for nm in ["H2-H2", "H2-He", "H2-H", "He-H"]:
            w = np.array(open(DATA + '/cia/' + nm + '_reform_11.txt').read().split(), float)
            n = int(w[0]); assert int(w[1]) == NB
            T = w[2:2+n]; kk = w[2+n+NB:2+n+NB+n*NB].reshape(n, NB)
            self.cia.append((T, kk))
        self.ray = []
        for nm in ["H2", "He", "H", "e-"]:
            with open(DATA + '/ray/Ray_' + nm + '_11.txt') as fh:
                fh.readline()
                self.ray.append(np.array(fh.read().split()[:NB], float))
        self.ray = np.array(self.ray)
        # stellar band fractions: blackbody at star_teff, tails folded (read_ck_continuum)
        fb = planck_below(self.wl*star_teff)
        s = fb[:-1] - fb[1:]
        s[0] += 1.0 - fb[0]
        s[-1] += fb[-1]
        if star_teff <= 0.0:
            # wasp121_0925: ck_star_teff <= 0 -> the band-flux file (ASCENDING wl, reversed
            # and normalised, exactly as read_ck_continuum does)
            v = np.loadtxt(_os.environ.get('RCE_SWFILE', DATA + '/sw_flux/sw_band_flux_W121_11.txt'))
            s = v[::-1]/v.sum()
        self.swf = s
        self.lam = 2.0/(1.0/self.wl[:-1] + 1.0/self.wl[1:])

    def continuum(self, T, pbar, rho):
        """ck_continuum, vectorised: (N,) -> (N, NB) cm^2/g."""
        T = np.asarray(T, float); pbar = np.asarray(pbar, float); rho = np.asarray(rho, float)
        iT, fT = tp_index(self.ceT, np.log10(T))
        iP, fP = tp_index(self.ceP, np.log10(pbar))
        c = self.ce
        vmr = ((1-fT)[:, None]*((1-fP)[:, None]*c[iT, iP] + fP[:, None]*c[iT, iP+1])
               + fT[:, None]*((1-fP)[:, None]*c[iT+1, iP] + fP[:, None]*c[iT+1, iP+1]))
        ntot = pbar*1e6/(KB*T)
        irho = 1.0/rho
        kc = np.zeros(T.shape + (NB,))
        i1 = [1, 1, 1, 2]; i2 = [1, 2, 3, 3]
        for s, (cT, ck) in enumerate(self.cia):
            n = len(cT)
            it = np.clip(np.searchsorted(cT, T, side='right') - 1, 0, n-2)
            ft = (T - cT[it])/(cT[it+1] - cT[it])
            lo = T <= cT[0]; hi = T >= cT[-1]
            it = np.where(lo, 0, np.where(hi, n-2, it)); ft = np.where(lo, 0., np.where(hi, 1., ft))
            nn = vmr[:, i1[s]]*ntot*vmr[:, i2[s]]*ntot*irho
            kc += ((1-ft)[:, None]*ck[it] + ft[:, None]*ck[it+1])*nn[:, None]
        for s in range(4):
            kc += self.ray[s][None, :]*(vmr[:, s+1]*ntot*irho)[:, None]
        # H- bf + ff (John 1988)
        lam0 = 1.6419
        Cbf = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]
        ff1 = np.array([[518.1021, 472.2636, -482.2089, 115.5291, 0.0, 0.0],
                        [-734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0],
                        [1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0],
                        [-479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0],
                        [93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0],
                        [-6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0]])
        ff2 = np.array([[0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830],
                        [0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170],
                        [0.0, -2054.2910, 8746.5230, -13651.1050, 8642.9700, -1863.8640],
                        [0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880],
                        [0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880],
                        [0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850]])
        T5040 = 5040.0/T
        nHm = vmr[:, 5]*ntot
        PenH = vmr[:, 4]*ntot*vmr[:, 3]*ntot*KB*T
        for b in range(NB):
            lam = self.lam[b]
            xbf = 0.0
            if lam < lam0:
                dk = 1.0/lam - 1.0/lam0
                sdk = np.sqrt(dk)
                fbf = sum(Cbf[n]*sdk**n for n in range(6))
                xbf = 1e-18*lam**3*(dk*sdk)*fbf
            sff = np.zeros_like(T)
            if lam >= 0.3645 or (0.1823 < lam < 0.3645):
                A = ff2 if lam >= 0.3645 else ff1
                for n in range(6):
                    tp = T5040**((n+2)/2.0)
                    sff += tp*(lam*lam*A[0, n] + A[1, n] + A[2, n]/lam + A[3, n]/lam**2
                               + A[4, n]/lam**3 + A[5, n]/lam**4)
            kc[:, b] += (xbf*nHm + 1e-29*sff*PenH)*irho
        return kc

    def kappa_line(self, T, pbar):
        """ck_kappa for all (b,g): (N,) -> (N, NB, NG)."""
        iT, fT = tp_index(self.lT, np.log10(T))
        iP, fP = tp_index(self.lP, np.log10(pbar))
        lk = self.lk
        l = ((1-fT)*((1-fP)*lk[:, :, iT, iP] + fP*lk[:, :, iT, iP+1])
             + fT*((1-fP)*lk[:, :, iT+1, iP] + fP*lk[:, :, iT+1, iP+1]))
        return np.moveaxis(10.0**l, -1, 0)

    def kr(self, T, pbar, rho, chunk=50000):
        """(kappa rho) per chain c = b*NG + g: (N,) -> (N, 88) float64, cm^-1."""
        T = np.ravel(T); pbar = np.ravel(pbar); rho = np.ravel(rho)
        out = np.zeros((T.size, NB*NG))
        for s in range(0, T.size, chunk):
            sl = slice(s, s+chunk)
            k = self.kappa_line(T[sl], pbar[sl]) + self.continuum(T[sl], pbar[sl], rho[sl])[:, :, None]
            out[sl] = (k*rho[sl, None, None]).reshape(-1, NB*NG)
        return out

    def chain_weight(self):
        """ckswf(b) * ckgw(g) per chain, (88,): sums to 1."""
        return (self.swf[:, None]*self.gw[None, :]).ravel()

    def rosseland_table(self):
        """ck_build_rosseland_table: log10 kappa_R on (lT, lP[bar]) grid."""
        TT, PP = np.meshgrid(10**self.lT, 10**self.lP, indexing='ij')
        T = TT.ravel(); pb = PP.ravel()
        iT, fT = tp_index(self.ceT, np.log10(T)); iP, fP = tp_index(self.ceP, np.log10(pb))
        c = self.ce[..., 0]
        mu = (1-fT)*((1-fP)*c[iT, iP] + fP*c[iT, iP+1]) + fT*((1-fP)*c[iT+1, iP] + fP*c[iT+1, iP+1])
        rho = pb*1e6*mu*1.6726e-24/(KB*T)
        kc = self.continuum(T, pb, rho)
        kl = 10.0**np.moveaxis(self.lk.reshape(NB, NG, -1), -1, 0)   # (N,b,g) at nodes
        Tp, Tm = 1.01*T, 0.99*T
        fp = np.array([self.planck_frac(t) for t in Tp]); fm = np.array([self.planck_frac(t) for t in Tm])
        wb = 5.670374419e-5*(Tp[:, None]**4*fp - Tm[:, None]**4*fm)
        inv = (self.gw[None, None, :]/(kl + kc[:, :, None])).sum(-1)
        wbp = np.where(wb > 0, wb, 0.0)
        tab = np.log10(wbp.sum(1)/(wbp*inv).sum(1))
        return tab.reshape(len(self.lT), len(self.lP))

    def planck_frac(self, T):
        fb = planck_below(self.wl*T)
        s = fb[:-1] - fb[1:]
        s[0] += 1.0 - fb[0]; s[-1] += fb[-1]
        return s


