"""premix k at the 18 lowP nodes + 1e-8 bar for every table T <= 6100 K, two SiO variants.
A: FastChem SiO; B: SiO x f(T) (table-implied, from fit3). FeII de-floored. Gauss-point g."""
import glob
import numpy as np
from multiprocessing import Pool
import premix_lib as p
fc = dict(np.load('/orion/ptmp/jinma/rsg_wind_1008/lowP/chem/fc_eq.npz'))
TAB = [s[0] for s in p.SP[:32]]
IPS = list(range(19))                      # 0..17 new nodes, 18 = 1e-8
assert abs(fc['P'][18]/1e-8-1) < 1e-9
LF = {}
for f in glob.glob('/orion/ptmp/jinma/rsg_wind_1008/dace/out/fit3_*.npy'):
    for T, v in np.load(f, allow_pickle=True)[0].items():
        LF[T] = v['lf']
TN = p.KT['T'][p.KT['T'] <= 6100.]


def fsio(T):
    lf = LF[T] if T <= 3100 else LF[3100.]
    return 0. if lf < -50 else 10**lf


def one(T):
    it = int(np.argmin(abs(fc['T']-T)))
    assert fc['T'][it] == T
    A = {}
    for s in TAB:
        a = p.load(s, T)[0].astype(np.float64)
        A[s] = p.defloor(a) if s == 'FeII' else a
    KA = np.zeros((19, p.NB, 8))
    KB = np.zeros((19, p.NB, 8))
    fs = fsio(T)
    for ip in IPS:
        X, mu = p.vmr(fc, it, ip, TAB)
        m = np.zeros(p.NMAX)
        for s in TAB:
            if s != 'SiO':
                m += X[s]*p.SPD[s][8]/mu*A[s]
        sio = X['SiO']*p.SPD['SiO'][8]/mu*A['SiO']
        KA[ip] = p.sortk(m+sio)[0]
        KB[ip] = p.sortk(m+fs*sio)[0]
    print(T, flush=True)
    return KA, KB, fs


if __name__ == '__main__':
    with Pool(12) as pool:
        out = pool.map(one, TN)
    np.savez('/orion/ptmp/jinma/rsg_wind_1008/dace/out/premix_lowP.npz', T=TN, P=fc['P'][:19],
             KA=np.array([o[0] for o in out]), KB=np.array([o[1] for o in out]),
             fsio=np.array([o[2] for o in out]))
    print('ALL')
