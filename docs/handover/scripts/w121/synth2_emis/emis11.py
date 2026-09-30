"""emis11.py: the new angle-dependent RT run with the SIMULATION's own 11-band Exo-FMS ck
opacity (hiT2 tables, the run's CE table, CIA, Rayleigh-as-absorption, H- bf/ff as the model).

  python emis11.py w1x|w10x
Output out/emis11_<arm>.npz:
  Fcol   (ncol, 11) top-face upward flux, 8-pt Gauss in mu, r^2 layer weighting (main)
  Fcol2  (ncol, 11) same with the 2-pt Gauss pair (the GCM's ck_nquad = 2)
  Fpp    (ncol, 11) 8-pt, plane-parallel (no r^2 weighting)
  Pang   (nph, 11)  disk flux x d^2, angle-dependent intensity
  Plam   (nph, 11)  disk flux x d^2, Lambertian from Fcol (I = F/pi), = the old projection
Bands in TABLE order (descending wavelength, as rt_Fb / olr files).
"""
import sys
import numpy as np
import rtcore as rc
import ck11

CK = {'w1x': dict(ktab='/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck/ck/Premixed_1x_g8_11_hiT2.txt',
                  cetab='/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck/CE_tables/FastChem_ck_1x_int_hiT2.txt',
                  ddir='/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck'),
      'w10x': dict(ktab='/viper/ptmp2/jinma/wasp121_0925/ckdata10/ck/Premixed_10x_g8_11_hiT2.txt',
                   cetab='/viper/ptmp2/jinma/wasp121_0925/ckdata10/CE_tables/FastChem_ck_10x_int_hiT2.txt',
                   ddir='/viper/ptmp2/jinma/wasp121_0925/ckdata10')}

if __name__ == '__main__':
    arm = sys.argv[1]
    st = np.load(rc.SD + 'state/%s.npz' % arm)
    T, p, rho, rf = st['T'], st['p'], st['rho'], st['rf']
    X, dOm = st['X'], st['dOm']
    ncol, nlay = T.shape
    ck = ck11.CK(ktab=CK[arm]['ktab'], cetab=CK[arm]['cetab'], ddir=CK[arm]['ddir'])
    kr = ck.kr(T.ravel(), p.ravel()/1e6, rho.ravel()).reshape(ncol, nlay, 11, 8)
    bp = rc.BandPlanck(ck.wl)
    B = bp(T)[..., ::-1]                        # (ncol, nlay, 11) table (descending) order
    dz = np.diff(rf)
    rcen = 0.5*(rf[1:] + rf[:-1])
    rtop = rf[-1]
    gx8, gw8 = rc.gauss01(8)
    gx2, gw2 = rc.gauss01(2)
    phi, obs = rc.phases(72)
    mup = X @ obs.T                               # (ncol, nph)
    Fcol = np.zeros((ncol, 11))
    Fcol2 = np.zeros((ncol, 11))
    Fpp = np.zeros((ncol, 11))
    Pang = np.zeros((len(phi), 11))
    gw = ck.gw
    for q in range(ncol):
        vis = mup[q] > 0
        mus = np.concatenate([gx8, gx2, mup[q, vis]])
        dtau = kr[q]*dz[:, None, None]
        S = (B[q]*rcen[:, None]**2)[:, :, None]
        J = rc.formal(dtau, rc.face_source(S, dtau), mus)          # (11, 8, nmu)
        Jg = np.einsum('bgm,g->bm', J, gw)
        Fcol[q] = 2*np.pi*(Jg[:, :8]*(gx8*gw8)).sum(1)/rtop**2
        Fcol2[q] = 2*np.pi*(Jg[:, 8:10]*(gx2*gw2)).sum(1)/rtop**2
        Pang[vis] += (Jg[:, 10:]*(mup[q, vis]*dOm[q])).T
        Spp = B[q][:, :, None]
        Jp = rc.formal(dtau, rc.face_source(Spp, dtau), gx8)
        Fpp[q] = 2*np.pi*(np.einsum('bgm,g->bm', Jp, gw)*(gx8*gw8)).sum(1)
    Plam = np.einsum('qb,qp->pb', Fcol*rtop**2/np.pi, np.maximum(mup, 0)*dOm[:, None])
    np.savez(rc.SD + 'out/emis11_%s.npz' % arm, Fcol=Fcol, Fcol2=Fcol2, Fpp=Fpp, Pang=Pang,
             Plam=Plam, phi=phi, edges=ck.wl, gid=st['gid'], k=st['k'], j=st['j'], rtop=rtop,
             dOm=dOm)
    print(arm, 'L_IR (8pt, r2) = %.4e  (2pt) %.4e  (pp) %.4e erg/s' % (
        (Fcol.sum(1)*dOm).sum()*rtop**2, (Fcol2.sum(1)*dOm).sum()*rtop**2,
        (Fpp.sum(1)*dOm).sum()*rtop**2))
