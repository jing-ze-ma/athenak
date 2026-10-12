"""fluxA.py prof.npz: per dump, luminosity decomposition (units L0 = L_in) at selected r/R_ph:
Lrad_co (comoving M1), Lconv = resolved fluctuation enthalpy (gas + rad) flux, Lkin, L_MLT = w F_MLT(IC) 4 pi r^2,
Lsum; plus band <v_r> [km/s], KE, photosphere (per-column tau=2/3 mean / max) in R_ph."""
import os, sys, numpy as np
RPH, RS = 2.70140777e13, 6.957e10
z = np.load(sys.argv[1]); r = z['r']; x = r/RPH
ic = np.loadtxt(os.environ['ICFILE'], comments='#')
f = ic[:, 6]; Fm = ic[:, 3]*f/np.maximum(1 - f, 1e-30)
LM = np.interp(r, ic[:, 0], 4*np.pi*ic[:, 0]**2*Fm)/3.411709e39
RS_, RT_ = float(os.environ.get('RAMP_S', '4.032e6')), float(os.environ.get('RAMP_T', '2.016e7'))
def w(t, s=RS_, T=RT_):
    return 1.0 if t < s else (0.0 if t >= s + T else 0.5*(1 + np.cos(np.pi*(t - s)/T)))
xs = [0.40, 0.45, 0.55, 0.63, 0.7, 0.8, 0.9, 0.95, 1.0]
ii = [int(np.argmin(abs(x - a))) for a in xs]
bands = [(0.34, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.03), (1.03, 1.5), (1.5, 3.0)]
every = int(sys.argv[2]) if len(sys.argv) > 2 else 1
print('r/R_ph =', ' '.join('%.3f' % x[i] for i in ii), '| L_MLT,IC/L0 there:', ' '.join('%.2f' % LM[i] for i in ii))
dv = np.diff(np.concatenate([[r[0] - z['dr'][0]/2], r + z['dr']/2])**3)/3   # per sr volume
for k in range(0, len(z['time']), every):
    t = float(np.ravel(z['time'])[k]); ww = w(t)
    Lr, Lc, Lk = z['Lrad_co'][k], z['Lconv_g'][k] + z['Lconv_r'][k], z['Lkin'][k]
    Lm = ww*LM; S = Lr + Lc + Lk + Lm
    vr = z['vr'][k]/1e5
    bv = ' '.join('%+6.1f' % vr[(x >= a) & (x < b)][np.argmax(abs(vr[(x >= a) & (x < b)]))] for a, b in bands)
    ke = (z['ke'][k]*dv)[x < 1.0].sum()*np.pi/3
    rp = z['rph_col'][k]*RS/RPH
    print('t %.3e w %.3f | KE(<R) %.2e | rph %.3f/%.3f | vr bands %s' % (t, ww, ke, rp[0], rp[2], bv))
    for nm, a in (('rad', Lr), ('conv', Lc), ('kin', Lk), ('mlt', Lm), ('sum', S)):
        print('   %-4s ' % nm + ' '.join('%6.3f' % a[i] for i in ii))
