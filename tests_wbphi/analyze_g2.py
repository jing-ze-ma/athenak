"""G2 analysis: max |v1| and the density drift of the three arms."""
import glob
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc                                      # noqa: E402

VMLT = 1.86e4
GAM1 = 1.66666


def load(path):
    d = bc.read_binary(path)
    vn = d['var_names']
    mb = d['mb_data']
    zs = []
    for n in range(d['n_mbs']):
        g = d['mb_geometry'][n]
        nx = d['nx1_out_mb']
        dz = (g[1] - g[0])/nx
        zs.append(g[0] + dz*(np.arange(nx) + 0.5))
    x1 = np.concatenate(zs)
    o = np.argsort(x1)
    out = {'t': d['time'], 'z': x1[o]}
    for k in vn:
        out[k] = np.concatenate([mb[k][n].ravel() for n in range(d['n_mbs'])])[o]
    return out


def arm(tag):
    fs = sorted(glob.glob('g2/%s/bin/*.bin' % tag))
    d0 = load(fs[0])
    cs = np.sqrt(GAM1*(GAM1 - 1.0)*d0['eint']/d0['dens'])
    vmax, tmax, zmax, csmax = 0.0, 0.0, 0.0, 1.0
    for f in fs:
        d = load(f)
        i = np.argmax(np.abs(d['velx']))
        if abs(d['velx'][i]) > vmax:
            vmax = abs(d['velx'][i])
            tmax = d['t']
            zmax = d['z'][i]
            csmax = cs[i]
    dl = load(fs[-1])
    ie = np.argmax(np.abs(dl['velx']))
    drift = np.abs(dl['dens']/d0['dens'] - 1.0)
    print('--- arm %s  (%d dumps, t_end = %.1f s)' % (tag, len(fs), dl['t']))
    print('  max|v1| over the run  = %.4e cm/s  at t = %.1f s, z = %.4e'
          % (vmax, tmax, zmax))
    print('     /c_s = %.3e   /v_MLT = %.3e' % (vmax/csmax, vmax/VMLT))
    print('  max|v1| at the end    = %.4e cm/s  at z = %.4e  /c_s = %.3e  /v_MLT = %.3e'
          % (abs(dl['velx'][ie]), dl['z'][ie],
             abs(dl['velx'][ie])/cs[ie], abs(dl['velx'][ie])/VMLT))
    print('  density drift: max %.4e  median %.4e  (at z = %.4e)'
          % (drift.max(), np.median(drift), dl['z'][np.argmax(drift)]))


for t in sys.argv[1:] or ['A', 'B', 'C']:
    arm(t)
