import glob, json, math, os, sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/rwdiag_1112/rw_amp1e6')
import run_radwave as rw
B = '/viper/ptmp2/jinma/rwdiag_1112/rw_wbs'
arms = ['FINAL', 'FIN', 'FINmslbe', 'FINmfct', 'FINmfmm', 'FINmwbs', 'FINmfold', 'FINpzold', 'FINpcrst', 'FINnp30']
for tau in (0.01, 0.1):
    print('tau', tau, '| arm | L1E 64 128 256 | pE 128->256 | Enm256 | L1rho256 | prho')
    for a in arms:
        G = {}
        for p in glob.glob(B + '/runs_%s/*/results.json' % a):
            r = json.load(open(p)); c = r['case']
            if c['tau'] != tau: continue
            bg, refs = rw.references(c); w, vec = refs['edd']
            x = np.array(r['x1']); e = np.array(r['e_end']); rho = np.array(r['rho_end']); t = r['t_end']
            A = float(c['extra']['amp']); k = bg.k
            eex = bg.e0*(1 + A*np.real(vec[3]/bg.e0*np.exp(1j*(k*x - w*t))))
            rex = bg.rho0*(1 + A*np.real(np.exp(1j*(k*x - w*t))))
            nE = A*abs(vec[3])
            G[c['N']] = (np.mean(np.abs(e - eex))/nE, np.mean(np.abs((e - e.mean()) - (eex - eex.mean())))/nE,
                         np.mean(np.abs(rho - rex))/(A*bg.rho0))
        if len(G) < 3: print(a, 'missing'); continue
        pE = math.log(G[128][0]/G[256][0])/math.log(2); pr = math.log(G[128][2]/G[256][2])/math.log(2)
        print('%-9s %.3e %.3e %.3e  %.2f  %.3e  %.3e  %.2f' % (a, G[64][0], G[128][0], G[256][0], pE, G[256][1], G[256][2], pr))
