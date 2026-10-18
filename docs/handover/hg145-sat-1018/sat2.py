"""sat2.py: as sat.py, but split the M1 lab-frame flux: L_adv = Omega r^2 [<4/3 E v_r> - 4/3 <E><v_r>] (radiation enthalpy carried by
the convective flow, part of m1_f1) and L_diff = Omega r^2 <F1> - L_adv. Ratios (L_conv,res + L_adv)/(L_in - L_diff) and (L_conv,res+L_adv)/L_in
per band; solid-angle weighted shell means as reduce.py."""
import os, sys, numpy as np
sys.path.insert(0, os.environ['VISPY']); import bin_convert as bc
RED = '/work/hdd/bivj/jma20/hg145_prod_1018/page/data'; BIN = '/work/hdd/bivj/jma20/hg145_prod_1018/bin'
HST = '/work/nvme/bivj/jma20/hg145_prod_1018/run/hegiant.user.hst'
RS = 6.957e10; OM = np.pi / 3
u = np.loadtxt(HST); _, keep = np.unique(u[::-1, 0], return_index=True); u = u[::-1][keep]
th, Lin = u[:, 0], u[:, 6]
LB = [(20, 30), (30, 60), (60, 100)]
print('# t[d] | L_adv/L_in | (L_conv,res+L_adv)/L_in | L_diff/L_in | (L_conv,res+L_adv)/(L_in-L_diff)   each at ' + ' '.join('%d-%d' % b for b in LB))
for i in sorted(set(list(range(0, 113, 8)) + [110, 111])):
    d = np.load('%s/red_%05d.npz' % (RED, i)); t = float(d['t']); L_in = Lin[np.argmin(abs(th - t))]
    h = bc.read_binary_as_athdf('%s/hegiant.hydro_w.%05d.bin' % (BIN, i)); m = bc.read_binary_as_athdf('%s/hegiant.m1.%05d.bin' % (BIN, i))
    vr = np.asarray(h['velx'], float); E = np.asarray(m['m1_e'], float); F1 = np.asarray(m['m1_f1'], float)
    thc = np.asarray(h['x2v'], float); r = np.asarray(h['x1v'], float)
    w = np.sin(thc)[None, :, None] * np.ones((vr.shape[0], 1, 1)); W = w.sum(axis=(0, 1))
    sm = lambda q: (q * w).sum(axis=(0, 1)) / W
    Ladv = OM * r ** 2 * (4. / 3.) * (sm(E * vr) - sm(E) * sm(vr)); Lrad = OM * r ** 2 * sm(F1); Ldiff = Lrad - Ladv
    Lc = OM * r ** 2 * (d['Fenth'] + d['Fke']); rr = r / RS
    cols = [[], [], [], []]
    for a, b in LB:
        s = (rr >= a) & (rr < b); A_, C_, D_ = Ladv[s].mean(), (Lc[s] + Ladv[s]).mean(), Ldiff[s].mean()
        cols[0].append(A_ / L_in); cols[1].append(C_ / L_in); cols[2].append(D_ / L_in); cols[3].append(C_ / (L_in - D_))
    print('%5.2f | ' % (t / 86400) + ' | '.join(' '.join('%+.3f' % x for x in c) for c in cols))
