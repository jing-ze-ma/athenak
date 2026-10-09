#!/usr/bin/env python3
"""Old vs new kappa_R, kappa_P along AG Car A/B columns: ICs (geos/icA, icB) and the evolved
1-D general-EOS columns (geos/col/colA, colB last dump, horizontal mean, T from eint by the
geos colprof Newton on the EOS dump).  Lookup = the code's RosselandTable (bilinear in log T,
log rho, clamped)."""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/planck_tools')
sys.path.insert(0, '/viper/ptmp2/jinma/lbv_1008/agcar/geos/scripts')
sys.dont_write_bytecode = True
from compare import read_repo  # noqa: E402

RS = 6.957e10
B = '/viper/ptmp2/jinma/lbv_1008/agcar'
TAB = {'old': (B + '/tables/rosseland_tops_gs98_x0.36_z0.02.txt',
               B + '/tables/planck_tops_gs98_x0.36_z0.02.txt'),
       'new': (B + '/tables_ext/rosseland_ext_gs98_x0.36_z0.02.txt',
               B + '/tables_ext/planck_ext_gs98_x0.36_z0.02.txt')}


def code_extend(tab, ld_ext=-21.0):
    """what he_star_m1's HsReadOpacityTable does with problem/he_opac_logd_min = -21 (the
    production key): extend to lower rho, log kappa linear in log rho with the slope of the
    table's first density interval clamped to [0, 1].  A no-op for the new table (starts at -21)."""
    lT, lD, K = tab
    dld = lD[1] - lD[0]
    if not ld_ext < lD[0] - 1e-9*dld:
        return tab
    nadd = int(np.ceil((lD[0] - ld_ext)/dld - 1e-9))
    sl = np.clip((K[:, 1] - K[:, 0])/dld, 0.0, 1.0)
    ext = K[:, :1] - sl[:, None]*dld*np.arange(nadd, 0, -1)[None, :]
    return lT, np.concatenate([lD[0] - dld*np.arange(nadd, 0, -1), lD]), np.hstack([ext, K])


T_ = {k: [code_extend(read_repo(f)) for f in v] for k, v in TAB.items()}


def look(tab, T, rho):
    lT, lD, K = tab
    x = np.clip(np.log10(T), lT[0], lT[-1])
    y = np.clip(np.log10(rho), lD[0], lD[-1])
    i = np.clip(np.searchsorted(lT, x, 'right') - 1, 0, len(lT) - 2)
    j = np.clip(np.searchsorted(lD, y, 'right') - 1, 0, len(lD) - 2)
    fx = (x - lT[i])/(lT[i+1] - lT[i])
    fy = (y - lD[j])/(lD[j+1] - lD[j])
    return 10**((1-fx)*(1-fy)*K[i, j] + fx*(1-fy)*K[i+1, j] + (1-fx)*fy*K[i, j+1]
                + fx*fy*K[i+1, j+1])


def cols():
    out = {}
    for X, rph in (('A', 388.3007), ('B', 101.2962)):
        ic = np.loadtxt(B + '/geos/ic%s/ic_agcar_%s_ge.txt' % (X, X))
        out['IC ' + X] = (ic[:, 0]/RS/rph, ic[:, 1], ic[:, 5], rph)
        D = B + '/geos/col/col%s' % X
        hb = sorted(glob.glob(D + '/bin/*.hydro_w.*.bin'))
        if hb:
            import bin_convert_agcar as bc
            import colprof_agcar_ge as _  # noqa: F401  (not imported: argv-driven)
    return out


def evolved(X, rph):
    import bin_convert_agcar as bc
    import make_ic_mlt_star_ge as mk
    EOS = mk.TableEOS(B + '/geos/eos/eos_dump.txt')
    D = B + '/geos/col/col%s' % X
    hb = sorted(glob.glob(D + '/bin/*.hydro_w.*.bin'))[-1]
    gf = glob.glob(D + '/*.x1grid.txt')[0]
    with open(gf) as fh:
        is_ = int(fh.readline().split()[3])
    grid = np.array([[float(v) for v in ln.split()] for ln in open(gf)
                     if not ln.startswith('#') and len(ln.split()) == 4])
    f = bc.read_binary(hb)
    mp = lambda v: np.mean([f['mb_data'][v][m].mean(axis=(0, 1))  # noqa: E731
                            for m in range(f['n_mbs'])], axis=0)
    rho, eint = mp('dens'), mp('eint')
    r = grid[is_:is_ + len(rho), 2]
    ic = np.loadtxt(B + '/geos/ic%s/ic_agcar_%s_ge.txt' % (X, X))
    lT = np.log10(np.exp(np.interp(r, ic[:, 0], np.log(ic[:, 5]))))
    lr = np.log10(np.maximum(rho, 1e-21))
    for _ in range(40):
        lT = np.clip(lT - (EOS.le(lT, lr) - np.log10(eint/rho))/np.maximum(
            EOS.le(lT, lr, dx=1), 1e-3), 3.6, 6.6)
    return r/RS/rph, rho, 10**lT, f['time']


def main():
    od = sys.argv[1]
    data = {}
    for X, rph in (('A', 388.3007), ('B', 101.2962)):
        ic = np.loadtxt(B + '/geos/ic%s/ic_agcar_%s_ge.txt' % (X, X))
        data['IC %s' % X] = (ic[:, 0]/RS/rph, ic[:, 1], ic[:, 5])
        try:
            x, rho, T, t = evolved(X, rph)
            data['col%s t=%.2e s' % (X, t)] = (x, rho, T)
        except Exception as e:  # noqa: BLE001
            print('evolved', X, 'skipped:', e)
    fig, ax = plt.subplots(len(data), 2, figsize=(12, 3.2*len(data)), squeeze=False)
    for k, (nm, (x, rho, T)) in enumerate(data.items()):
        kRo, kPo = look(T_['old'][0], T, rho), look(T_['old'][1], T, rho)
        kRn, kPn = look(T_['new'][0], T, rho), look(T_['new'][1], T, rho)
        print('== %s: x=r/R_ph %.3f..%.3f' % (nm, x.min(), x.max()))
        for a, b in ((0.0, 0.9), (0.9, 1.0), (1.0, 1.1), (1.1, 1.5), (1.5, 2.0), (2.0, 3.2)):
            m = (x >= a) & (x < b)
            if not m.any():
                continue
            dR = np.log10(kRn[m]/kRo[m])
            dP = np.log10(kPn[m]/kPo[m])
            print('  r/R_ph %.1f-%.1f n %5d  T %6.0f-%6.0f K rho %.1e-%.1e | '
                  'kR old %.3g-%.3g new %.3g-%.3g  dlog max %+.2f/%+.2f | '
                  'kP old %.3g-%.3g new %.3g-%.3g dlog %+.2f/%+.2f | bitwise R %d/%d' % (
                      a, b, m.sum(), T[m].min(), T[m].max(), rho[m].min(), rho[m].max(),
                      kRo[m].min(), kRo[m].max(), kRn[m].min(), kRn[m].max(), dR.min(), dR.max(),
                      kPo[m].min(), kPo[m].max(), kPn[m].min(), kPn[m].max(), dP.min(), dP.max(),
                      (kRn[m] == kRo[m]).sum(), m.sum()))
        for c, (kn, ko, lab) in enumerate(((kRn, kRo, 'kappa_R'), (kPn, kPo, 'kappa_P'))):
            a = ax[k, c]
            a.loglog(x, ko, label='old (TOPS + code rho-extension, he_opac_logd_min -21)', lw=1.2)
            a.loglog(x, kn, '--', label='new (ext)', lw=1.2)
            if c == 0:
                a.axhline(0.2677, color='gray', lw=0.6, ls=':')
            a.axvspan(0.9, 3.0, color='0.9')
            a.set_title('%s  %s' % (nm, lab), fontsize=9)
            a.set_xlabel('r / R_ph')
            a2 = a.twinx()
            a2.semilogy(x, T, color='C3', lw=0.6, alpha=0.6)
            a2.set_ylabel('T [K] (red)', fontsize=7)
            a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(od + '/kappa_columns.png', dpi=110)
    print('wrote', od + '/kappa_columns.png')


if __name__ == '__main__':
    main()
