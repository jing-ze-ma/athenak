# usage: ring.py <run dir> [t0] [t1]: radial-ringing period (FFT of <v_r>) and envelope growth rates
import sys, glob, numpy as np
d = sys.argv[1]; gam = 5/3
f = glob.glob(d + '/*.user.hst')[0]
hdr = [l for l in open(f) if l.startswith('#')][-1]
names = [s.split('=')[1] for s in hdr[1:].split()]
a = np.atleast_2d(np.loadtxt(f)); c = {n: a[:, i] for i, n in enumerate(names)}
t = c['time']; t0 = float(sys.argv[2]) if len(sys.argv) > 2 else t[0]; t1 = float(sys.argv[3]) if len(sys.argv) > 3 else t[-1]
m = (t >= t0) & (t <= t1); t = t[m]
cs = np.sqrt(gam*c['PV_int'][m]/c['M_int'][m]); vr = c['Mr_int'][m]/c['M_int'][m]/cs
KE = c['KE_int'][m]; KEr = c['KEr_int'][m]; KEl = np.maximum(KE - KEr, 1e-300)
PV = c['PV_int'][m]
tu = np.linspace(t[0], t[-1], len(t)); y = np.interp(tu, t, vr); y -= y.mean()
P = np.abs(np.fft.rfft(y*np.hanning(len(y))))**2; fr = np.fft.rfftfreq(len(y), tu[1]-tu[0])
k = np.argmax(P[1:]) + 1
print('t %.0f..%.0f  <v_r>/cs peak period %.1f s (freq %.3e Hz); 2nd peak %.1f s' % (t[0], t[-1], 1/fr[k], fr[k],
      1/fr[np.argsort(P[1:])[-2]+1]))
nb = 6; e = np.array_split(np.arange(len(t)), nb)
for g in e:
    print('  t %7.0f-%7.0f  Mach_tot %.3e  Mach_r %.3e  Mach_lat %.3e  |<vr>|max/cs %.2e' % (t[g[0]], t[g[-1]],
          np.sqrt(2*KE[g].mean()/(gam*PV[g].mean())), np.sqrt(2*KEr[g].mean()/(gam*PV[g].mean())),
          np.sqrt(2*KEl[g].mean()/(gam*PV[g].mean())), abs(vr[g]).max()))
tm = np.array([t[g].mean() for g in e])
for nm, arr in (('KE_r', KEr), ('KE_lat', KEl)):
    lk = np.log([arr[g].mean() for g in e]); p = np.polyfit(tm, lk, 1)
    print('  %s amplitude e-folding time %.3e s (growth rate of amplitude %.3e /s)' % (nm, 2/p[0], p[0]/2))
