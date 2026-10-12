"""atm.py RPH GM run... : IC-atmosphere drift metric (reatm-1011).
Per dump: shell-mean <v_r> [km/s] (max |.| over cells in 1.00-1.05 R_ph and per band), mass-weighted band <v_r>,
band mass / mass0 - 1, and a = <v_r>/t in units of g (band centre)."""
import sys, glob, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import bin_convert as bc
R = float(sys.argv[1]); GM = float(sys.argv[2])
bands = [(0.97, 1.0), (1.0, 1.01), (1.01, 1.02), (1.02, 1.03), (1.03, 1.05)]
for d in sys.argv[3:]:
    fs = sorted(glob.glob(d + '/bin/agcar3d.hydro_w.*.bin'))
    if not fs: continue
    h0 = bc.read_binary_as_athdf(fs[0]); r = np.asarray(h0['x1v'], float); xf = np.asarray(h0['x1f'], float)
    vol = (xf[1:]**3 - xf[:-1]**3)/3.0; xr = r/R
    rho0 = np.asarray(h0['dens'], float).mean(axis=(0, 1))
    print('#### %s  (bands %s; per band: <v_r>_m km/s, dM/M0, a/g)' % (d.rstrip('/').split('/')[-1], ' '.join('%.2f-%.2f' % b for b in bands)))
    w5 = (xr >= 1.0) & (xr < 1.05)
    for f in fs:
        h = bc.read_binary_as_athdf(f); c = int(h['NumCycles']); t = float(h['Time'])
        rho = np.asarray(h['dens'], float); vx = np.asarray(h['velx'], float)
        rm = rho.mean(axis=(0, 1)); vm = vx.mean(axis=(0, 1))/1e5; pm = (rho*vx).mean(axis=(0, 1))
        s = []
        for a, b in bands:
            w = (xr >= a) & (xr < b)
            M = np.sum(rm[w]*vol[w]); M0 = np.sum(rho0[w]*vol[w]); P = np.sum(pm[w]*vol[w])
            vb = P/M; g = GM/(R*0.5*(a + b))**2
            s.append('%+7.3f %+6.3f %+.1e' % (vb/1e5, M/M0 - 1, vb/t/g if t > 0 else 0.0))
        k = np.argmax(np.abs(vm[w5]))
        print('c%-5d t %.3e max|<vr>|(1-1.05) %+7.3f at %.4f | %s' % (c, t, vm[w5][k], xr[w5][k], ' | '.join(s)))
