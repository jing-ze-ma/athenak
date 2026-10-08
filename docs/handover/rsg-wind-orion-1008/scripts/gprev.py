"""thin Gamma_F golden16 on T 1000-3500 x rho 1e-17..1e-14 (lowP/preview/preview.py part a)."""
import sys, os
sys.path.insert(0, '/orion/ptmp/jinma/rsg_wind_1008/lowP/preview')
sys.path.insert(1, '/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/scripts')
import numpy as np  # noqa
import rsglib as L  # noqa
tag = sys.argv[1]
st = L.star('golden16')
Tg = np.arange(1000., 3501., 250.)
lrg = np.arange(-17., -13.99, 0.5)
G = np.zeros((len(Tg), len(lrg))); Pb = np.zeros_like(G); KL = np.zeros((len(Tg), len(lrg), L.NB))
for i, T in enumerate(Tg):
    for j, lr in enumerate(lrg):
        s = L.state(T, 10**lr, 'clamp')
        G[i, j] = L.kappa_thin_flux(s['K'], st['Teff'])/st['kE']
        Pb[i, j] = s['P']*1e-6
        KL[i, j] = (L.GW*s['kline']).sum(-1)
np.savez(f'/orion/ptmp/jinma/rsg_wind_1008/dace/out/gprev_{tag}.npz', Tg=Tg, lrg=lrg, G=G, Pb=Pb, KL=KL)
print(tag, 'done; Pbar range %.1e..%.1e' % (Pb.min(), Pb.max()))
