"""state.py: extract the rot-300 column state of w1x / w10x into a compact npz (read-only on the dumps).

  python state.py w1x|w10x
Output state/<arm>.npz: T, p, rho (ncol, nx1) cell values; rf (nx1+1) face radii, rc (nx1) centres;
X (ncol,3) column unit vector (dhjcs frame: substellar = -x); dOm (ncol) solid angle;
gid/k/j index of each column (column order = (m, k, j) of the bin).
"""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts/w121/synth_rot300')
import dhjcs   # noqa: E402
from synth import geometry   # noqa: E402

P = '/viper/ptmp2/jinma/w121prod_0929/'
ARMS = {
    'w1x': dict(bin=P + 'w1x/bin/dhj.hydro_w.00150.bin', eos=P + 'synth_rot300/diag/eos_table_1x.txt',
                R0=1.126575e10, R1=1.631037e10, N=76,
                C=[-0.379629, -1.797216, 7.873416, -79.711522, 347.885804, -753.943230,
                   784.221097, -314.370119]),
    'w10x': dict(bin=P + 'w10x/bin/dhj.hydro_w.00150.bin',
                 eos=P + 'synth_rot300_10x/diag/eos_table_10x.txt',
                 R0=1.116018e10, R1=1.530921e10, N=74,
                 C=[-0.090789, -4.624280, 43.684697, -282.659406, 932.446804, -1642.475346,
                    1459.195482, -514.767866]),
}


def stretch(xi, a):
    u = xi.copy()
    xk = xi.copy()
    for ck in a['C']:
        u += ck*xk*(1 - xi)
        xk = xk*xi
    return a['R0'] + (a['R1'] - a['R0'])*u


if __name__ == '__main__':
    arm = sys.argv[1]
    a = ARMS[arm]
    raw = dhjcs.bin_convert.read_binary(a['bin'])
    X, _, _, dOm = geometry(raw)
    mb = raw['mb_data']
    rho = np.asarray(mb['dens'], float)
    T, p = dhjcs.EOS(a['eos']).invert(rho, np.asarray(mb['eint'], float))
    nmb, nk, nj, ni = rho.shape
    assert ni == a['N']
    N = a['N']
    rf = stretch(np.arange(N + 1)/N, a)
    rc = stretch((np.arange(N) + 0.5)/N, a)
    assert np.all(np.diff(rf) > 0) and np.all((rc > rf[:-1]) & (rc < rf[1:]))
    M, K, J = np.meshgrid(np.arange(nmb), np.arange(nk), np.arange(nj), indexing='ij')
    np.savez(P + 'synth2_emis/state/%s.npz' % arm, T=T.reshape(-1, ni), p=p.reshape(-1, ni),
             rho=rho.reshape(-1, ni), rf=rf, rc=rc, X=X.reshape(-1, 3), dOm=dOm.ravel(),
             gid=M.ravel(), k=K.ravel(), j=J.ravel(), time=raw['time'])
    print(arm, 'time', raw['time'], 'rot', raw['time']/1.101535e5, 'T range', T.min(), T.max(),
          'p range [bar]', p.min()/1e6, p.max()/1e6, 'top p [bar] min/max',
          p[..., -1].min()/1e6, p[..., -1].max()/1e6, 'ncol', rho[..., 0].size)
