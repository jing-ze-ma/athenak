"""valid.py: validation 7a (static isothermal -> zero shift; solid-body rotation -> analytic).
Analytic rotation prediction: the static limb spectrum shifted by Omega b_eff cos(theta)
(b_eff = absorption-weighted impact parameter from the static run) and averaged over the
position angles of the limb; measured with the same Gaussian fit as ana.py (R = 0)."""
import numpy as np
from ana import gfit
W = '/viper/ptmp2/jinma/w121prod_0929/synth2_trans/out/'
OM = 5.704026e-05
mw = np.load(W + '1x_template.npz')['mw']
s0, s1 = np.load(W + '1x_iso0.npz'), np.load(W + '1x_isorot.npz')
th = np.deg2rad(np.arange(0, 360, 2.0) + 1.0)
for sp, U, key, EWk in (('Fe CCF', s0['UW'], 'fe', 0), ('Na D2', s0['UNA'], 'na', 1)):
    def spec(z):
        if key == 'fe':
            return np.einsum('w,elwu->lu', mw, z['Afe'] - z['Afe_c'][..., None])
        return (z['Ana'] - z['Ana_c'][..., None])[0, :, 0]
    st, ro = spec(s0), spec(s1)
    EW, bb = s0['CEb'][0, EWk], s0['bb'][0]   # line-centre absorption weights
    beff = (EW*bb).sum()/EW.sum()
    for limb, sg in ((0, 1), (1, -1)):
        tl = th[sg*np.cos(th) > 0]
        pred = np.mean([np.interp(U, U + OM*beff*np.cos(t), st[limb]) for t in tl], 0)
        print('%s %s  static %+6.3f   rotating: code %+6.3f  analytic %+6.3f  km/s'
              '  (b_eff %.4e cm, Omega b_eff (2/pi) %.3f)' % (
                  sp, ('morning', 'evening')[limb], gfit(U, st[limb])[0], gfit(U, ro[limb])[0],
                  gfit(U, pred)[0], beff, OM*beff*2/np.pi/1e5))
    predb = np.mean([np.interp(U, U + OM*beff*np.cos(t), st.sum(0)/2) for t in th], 0)*2
    print('%s both: code centre %+6.3f sigma %.3f ; analytic centre %+6.3f sigma %.3f' % (
        sp, *gfit(U, ro.sum(0))[:2], *gfit(U, predb)[:2]))
