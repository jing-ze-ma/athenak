import sys, glob, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert

def gridfrac(f):
    """fraction of kinetic energy in the odd-even (grid-scale) modes of the dump."""
    r = bin_convert.read_binary(f); mb = r['mb_data']
    rho = np.asarray(mb['dens'])
    ke = 0
    kg = 0
    for v in ('velx', 'vely', 'velz'):
        u = np.asarray(mb[v])
        ke += np.sum(rho*u*u)
        # project onto (-1)^(i+j+k)
        nm, nk, nj, ni = u.shape
        k, j, i = np.meshgrid(np.arange(nk), np.arange(nj), np.arange(ni), indexing='ij')
        sg = ((i + j + k) % 2 == 0)*2.0 - 1.0
        amp = np.sum(u*sg[None], axis=(1, 2, 3))/(nk*nj*ni)
        kg += np.sum(rho*(amp[:, None, None, None]*sg[None])**2)
    return r['time'], ke, kg/max(ke, 1e-300)

print('%-13s %8s %10s %10s %10s %10s   %s' % ('run', 't_end', 'KE(0)', 'KE(min)', 'KE(end)', 'KE(end)/KE(0)', 'grid-scale KE fraction at dumps'))
for d in sys.argv[1:]:
    h = np.atleast_2d(np.loadtxt(glob.glob(d + "/*.hst")[0]))
    ke = h[:, 7] + h[:, 8] + h[:, 9]
    fr = [gridfrac(f) for f in sorted(glob.glob(d + '/bin/*.bin'))]
    fin = np.isfinite(ke).all()
    print('%-13s %8.1f %10.3e %10.3e %10.3e %10.3g   %s %s' % (
        d, h[-1, 0], ke[0], ke.min(), ke[-1], ke[-1]/ke[0],
        ' '.join('%.2f' % g[2] for g in fr), '' if fin else 'NON-FINITE'))
