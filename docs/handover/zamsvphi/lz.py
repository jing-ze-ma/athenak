# ZAMSVPHI L_z reductions of M35ZAMS bins: where the thin layer's opposite L_z sits, and its pattern
# usage: python lz.py BINDIR [BASENAME] [READER_DIR]   (prints one block per dump)
import sys, glob, types, numpy as np
sys.modules.setdefault('h5py', types.ModuleType('h5py'))
D = sys.argv[1]; B = sys.argv[2] if len(sys.argv) > 2 else 'm35zams'
sys.path.insert(0, sys.argv[3] if len(sys.argv) > 3 else
                '/resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/vis/python')
import bin_convert as bc
C = 2.99792458e10; KES = 0.34
BANDS = (('core i<300', 0, 300), ('i300-319', 300, 320), ('thin i320-331', 320, 332),
         ('i332-343', 332, 344), ('floor i344+', 344, 100000))
NB = 8   # lateral bins for the theta / phi pattern of the thin layer
for f in sorted(glob.glob(D + '/' + B + '.vet.*.bin')):
    k = f.split('.')[-2]
    try: v = bc.read_binary(f); w = bc.read_binary(f.replace('.vet.', '.hydro_w.'))
    except Exception as e: print(k, 'read error', e); continue
    x1f = np.array(v['mb_x1f']); x2f = np.array(v['mb_x2f']); x3f = np.array(v['mb_x3f'])
    nb = len(v['mb_data']['vet_f1']); n1 = x1f.shape[1] - 1
    lzg = np.zeros(n1); lzr = np.zeros(n1); tq = np.zeros(n1); ms = np.zeros(n1)
    rey = np.zeros(n1 + 1); nfl = np.zeros(n1); mvt = np.zeros(n1)
    j0 = x2f[:, 0].min(); j1 = x2f[:, -1].max(); k0 = x3f[:, 0].min(); k1 = x3f[:, -1].max()
    pj = np.zeros(NB); pk = np.zeros(NB); pjs = np.zeros(NB); pks = np.zeros(NB)
    for b in range(nb):
        rf = x1f[b]; r = 0.5*(rf[1:] + rf[:-1]); th = 0.5*(x2f[b, 1:] + x2f[b, :-1])
        ph = 0.5*(x3f[b, 1:] + x3f[b, :-1])
        r3 = (rf[1:]**3 - rf[:-1]**3)/3.0; dc = np.cos(x2f[b, :-1]) - np.cos(x2f[b, 1:])
        dp = x3f[b, 1:] - x3f[b, :-1]
        V = dp[:, None, None]*dc[None, :, None]*r3[None, None, :]
        R = r[None, None, :]*np.sin(th)[None, :, None]
        g = lambda d, n: np.array(d['mb_data'][n][b], dtype=float)
        rho = g(w, 'dens'); vr = g(w, 'velx'); vt = g(w, 'vely'); vp = g(w, 'velz'); f3 = g(v, 'vet_f3')
        lg = rho*vp*R*V
        lzg += lg.sum(axis=(0, 1)); ms += (rho*V).sum(axis=(0, 1))
        lzr += (f3/C**2*R*V).sum(axis=(0, 1)); mvt += (rho*vt*V).sum(axis=(0, 1))
        tq += (KES*rho*f3/C*R*V).sum(axis=(0, 1))   # electron-scattering lower bound of the radiation torque
        nfl += (rho < 1.0e-14).sum(axis=(0, 1))
        # Reynolds L_z flux through the inner x1 face of each cell (cell-centre product, face mean)
        q = rho*vr*vp*R
        A = dp[:, None]*dc[None, :]
        rey[1:n1] += ((0.5*(q[:, :, 1:] + q[:, :, :-1]))*A[:, :, None]*(rf[None, None, 1:-1]**2)).sum(axis=(0, 1))
        sl = lg[:, :, 320:332].sum(axis=2); sa = np.abs(lg[:, :, 320:332]).sum(axis=2)
        ib = np.minimum(((th - j0)/(j1 - j0)*NB).astype(int), NB - 1)
        kb = np.minimum(((ph - k0)/(k1 - k0)*NB).astype(int), NB - 1)
        for jj in range(len(th)):
            pj[ib[jj]] += sl[:, jj].sum(); pjs[ib[jj]] += sa[:, jj].sum()
        for kk in range(len(ph)):
            pk[kb[kk]] += sl[kk, :].sum(); pks[kb[kk]] += sa[kk, :].sum()
    print('== dump %s t %.6g cycle %d' % (k, v['time'], v['cycle']))
    print('  band            M[g]        Lz_gas      Lz_rad      tq_es[dLz/dt]  <v_phi>_m   n(rho<1e-14)')
    for lab, a, z in BANDS:
        s = slice(a, min(z, n1)); Rm = 0.5*(x1f[0, a] + x1f[0, min(z, n1)])
        print('  %-14s %11.4e %11.4e %11.4e %11.4e %11.4e %8d' % (lab, ms[s].sum(), lzg[s].sum(),
              lzr[s].sum(), tq[s].sum(), lzg[s].sum()/max(ms[s].sum(), 1e-300)/Rm, nfl[s].sum()))
    print('  thin <v_theta>_m %.4e' % (mvt[320:332].sum()/ms[320:332].sum()))
    print('  TOTAL gas %.6e rad %.6e   (gas+rad %.6e)' % (lzg.sum(), lzr.sum(), lzg.sum() + lzr.sum()))
    print('  Reynolds L_z flux (outward, g cm^2/s^2) at x1 faces i=300 %.4e 320 %.4e 332 %.4e 344 %.4e'
          % (rey[300], rey[320], rey[332], rey[344]))
    print('  thin-layer Lz by theta bin (net/abs):', ' '.join('%+.2f' % (pj[i]/max(pjs[i], 1e-300)) for i in range(NB)),
          '| share of net:', ' '.join('%+.2f' % (pj[i]/pj.sum()) for i in range(NB)))
    print('  thin-layer Lz by phi bin   (net/abs):', ' '.join('%+.2f' % (pk[i]/max(pks[i], 1e-300)) for i in range(NB)),
          '| share of net:', ' '.join('%+.2f' % (pk[i]/pk.sum()) for i in range(NB)))
    ii = np.arange(290, min(n1, 384), 2)
    print('  per-shell (i, Lz_gas, Lz_rad, M):', ' '.join('%d:%.3e/%.2e/%.2e' % (i, lzg[i], lzr[i], ms[i]) for i in ii))
