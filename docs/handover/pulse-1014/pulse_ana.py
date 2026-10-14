"""pulse_ana.py MODE PATH [OUT]   (PU pulse-1014; python3 + numpy)
MODE run     : PATH = run dir with bin/*.hydro_w.*.bin  -> OUT (default PATH/pulse.npz) + table on stdout
MODE caltech : PATH = dir of e1div_prof.py prof_*.npz   -> same table/npz from their sin-weighted shell means
Per dump, from the sin(theta)-weighted shell means <rho>(r), <v_r>(r) (as e1div_prof.py):
  vrV  = r^2 dr-weighted mean of <v_r> over 1.0-1.2 R_ph (km/s)
  vrM  = <rho> r^2 dr-weighted mean of <v_r> over 1.0-1.2 R_ph (km/s)
  vr105, vr110 = <v_r> at 1.05, 1.10 R_ph (km/s);  d105, d110 = <rho>/<rho>(t=0) there
  KEe  = 0.5 <rho> <v_r>^2 r^2 dr integrated over 0.9-1.3 R_ph (shell-mean radial KE per sr, arbitrary units: compare runs)
  rsurf = outermost r/R_ph where <rho> > 1e-12 (the photosphere/inflated layer edge)
"""
import glob, os, sys
import numpy as np
RPH = 2.70140777e13
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def shell_from_bin(f):
    import bin_convert as bc
    h = bc.read_binary_as_athdf(f)
    th = np.asarray(h['x2v'], float); w = np.sin(th)[None, :, None]; w = w/np.sin(th).sum()/h['dens'].shape[0]
    rho = np.asarray(h['dens'], float); vr = np.asarray(h['velx'], float)
    vt = np.sqrt(np.asarray(h['vely'], float)**2 + np.asarray(h['velz'], float)**2)
    return (float(h['Time']), np.asarray(h['x1v'], float), (rho*w).sum(axis=(0, 1)), (vr*w).sum(axis=(0, 1)),
            (rho*vr*w).sum(axis=(0, 1)), (rho*vt**2*w).sum(axis=(0, 1)), (rho*vr**2*w).sum(axis=(0, 1)))


def metrics(t, r, rho, vr, rho0, rvr=None, rvt2=None, rvr2=None):
    x = r/RPH; dr = np.gradient(r); wv = r**2*dr
    s = (x >= 1.0) & (x <= 1.2)
    vrV = (vr[s]*wv[s]).sum()/wv[s].sum()/1e5
    mflux = (rvr if rvr is not None else rho*vr)
    vrM = (mflux[s]*wv[s]).sum()/(rho[s]*wv[s]).sum()/1e5
    o = dict(t=t, vrV=vrV, vrM=vrM)
    for xx, n in ((1.05, '105'), (1.10, '110')):
        i = np.argmin(abs(x - xx)); o['vr' + n] = vr[i]/1e5; o['d' + n] = rho[i]/rho0[i]
    e = (x >= 0.9) & (x <= 1.3)
    o['KEe'] = (0.5*(rvr2 if rvr2 is not None else rho*vr**2)[e]*wv[e]).sum()
    o['KEt'] = (0.5*rvt2[e]*wv[e]).sum() if rvt2 is not None else np.nan
    j = np.where(rho > 1e-12)[0]; o['rsurf'] = x[j.max()] if len(j) else np.nan
    return o


def main():
    mode, path = sys.argv[1], sys.argv[2]
    rows, prof = [], {}
    if mode == 'run':
        fs = sorted(glob.glob(os.path.join(path, 'bin', '*.hydro_w.*.bin')))
        rho0 = None
        for f in fs:
            t, r, rho, vr, rvr, rvt2, rvr2 = shell_from_bin(f)
            if rho0 is None: rho0 = rho
            rows.append(metrics(t, r, rho, vr, rho0, rvr, rvt2, rvr2)); prof.setdefault('vr', []).append(vr); prof.setdefault('rho', []).append(rho)
        out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(path, 'pulse.npz')
    else:
        fs = sorted(glob.glob(os.path.join(path, 'prof_*.npz')))
        rho0 = None
        for f in fs:
            d = np.load(f); r = d['r']; rho = d['rho_mean']; vr = d['vr_mean']
            if rho0 is None: rho0 = rho
            rows.append(metrics(float(d['time']), r, rho, vr, rho0)); prof.setdefault('vr', []).append(vr); prof.setdefault('rho', []).append(rho)
        out = sys.argv[3] if len(sys.argv) > 3 else 'caltech_pulse.npz'
    rows.sort(key=lambda o: o['t'])
    ks = ['t', 'vrV', 'vrM', 'vr105', 'vr110', 'd105', 'd110', 'KEe', 'KEt', 'rsurf']
    print('# ' + ' '.join('%10s' % k for k in ks))
    for o in rows: print('  ' + ' '.join('%10.3e' % o[k] for k in ks))
    np.savez(out, r=r, **{k: np.array([o[k] for o in rows]) for k in ks}, **{k: np.array(v) for k, v in prof.items()})


if __name__ == '__main__':
    main()
