import numpy as np, sys
TROT, LSW = 110154.0, None
def load(fn):
    with open(fn) as f:
        f.readline(); h = f.readline()
    names = [s.split('=')[1] for s in h.split()[1:]]
    d = np.loadtxt(fn); return {k: d[:, i] for i, k in enumerate(names)}
for run in sys.argv[1:]:
    u = load(f'{run}/dhj.user.hst'); m = load(f'{run}/dhj.hydro.hst')
    print(f'== {run}  (per window, units Lsw_abs mean; Lsw_abs = {u["Lsw_abs"][len(u["Lsw_abs"])//2]:.4e} erg/s)')
    for lo, hi in ((0, 1), (1, 2), (2, 2.6)):
        im = np.where((m['time'] >= lo*TROT) & (m['time'] <= hi*TROT))[0]
        iu = np.where((u['time'] > lo*TROT) & (u['time'] <= hi*TROT))[0]
        L = u['Lsw_abs'][iu].mean()
        dE = (m['tot-E'][im[-1]] - m['tot-E'][im[0]])/(m['time'][im[-1]] - m['time'][im[0]])
        mn = lambda c: u[c][iu].mean()
        res = (dE - mn('Etot_bot') + mn('Etot_top') - mn('Efloor'))/L
        print(f' rot {lo}-{hi}: dE/dt {dE/L:+.4e}  Lir_top {mn("Lir_top")/L:.5f}  Lsw_refl {mn("Lsw_refl")/L:.4f}  Etot_top {mn("Etot_top")/L:+.4e}'
              f'  Etot_bot {mn("Etot_bot")/L:+.4e}  Efloor {mn("Efloor")/L:+.3e}  Mdot_top {mn("Mdot_top"):+.3e} g/s  resid {res:+.3e}  mass drift {m["mass"][im[-1]]/m["mass"][0]-1:+.2e}')
