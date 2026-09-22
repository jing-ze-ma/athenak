"""Q2: T(p) day / terminator / night per arm, and the radial T profile in the
eta-capped shell.  Read-only."""
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert                       # noqa: E402
import dhjcs                             # noqa: E402
import common as C                       # noqa: E402

EOS = C.EOS2()
BARS = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0])
LEV = np.log10(BARS * 1e6)               # log10 p in barye


def regions(lon):
    a = np.abs(np.degrees(lon))
    return {'day': a < 30, 'term': (a > 75) & (a < 105), 'night': a > 150}


def profile(path):
    raw = bin_convert.read_binary(path)
    mb = raw['mb_data']
    rho = np.asarray(mb['dens'], dtype=np.float64)
    T, p, _ = EOS.invert(rho, np.asarray(mb['eint'], dtype=np.float64))
    geo = dhjcs.cs_geometry(raw)
    g = C.cs_metric(raw)
    lp = np.log10(p)
    Tl = dhjcs.level_interp(T, lp, LEV)          # (nlev, nmb, nk, nj)
    reg = regions(geo['lon'])
    out = {}
    for name, msk in reg.items():
        m3 = np.broadcast_to(msk, Tl.shape[1:])
        vals = np.where(m3[None], Tl, np.nan)
        out[name] = np.nanmean(vals.reshape(len(LEV), -1), axis=1)
    # horizontally averaged T(r) per region, for the eta-capped shell question
    out['r'] = g['rc'] / 9.44e9
    for name, msk in reg.items():
        out['Tr_' + name] = np.nanmean(T[msk], axis=0)
    out['t'] = raw['time']
    out['rho_r'] = np.nanmean(rho.reshape(-1, rho.shape[-1]), axis=0)
    return out


def main():
    prof = {}
    for name, d in C.ARMS1:
        f = sorted(glob.glob(d + '/bin/*.bin'))[-1]
        prof[name] = profile(f)

    fh = open(C.OUT + 'tp_tables.md', 'w')
    fh.write('# T(p) by region (tp.py)\n\n'
             'Horizontal mean of T interpolated to each pressure level; '
             'day = |lon_ss|<30 deg, terminator = 75-105 deg, night = |lon_ss|>150 deg. '
             'Each arm at ITS OWN last dump (times differ, see cfl_tables.md).\n')
    for reg in ['day', 'term', 'night']:
        fh.write('\n### %s: T [K]\n\n| p [bar] | ' % reg
                 + ' | '.join(n for n, _ in C.ARMS1)
                 + ' | max spread |\n|---|' + '---|' * (len(C.ARMS1) + 1) + '\n')
        for i, b in enumerate(BARS):
            v = [prof[n][reg][i] for n, _ in C.ARMS1]
            fh.write('| %.0e | ' % b + ' | '.join('%.0f' % x for x in v)
                     + ' | %.0f |\n' % (np.nanmax(v) - np.nanmin(v)))
    fh.write('\n### day - night contrast [K]\n\n| p [bar] | '
             + ' | '.join(n for n, _ in C.ARMS1) + ' |\n|---|'
             + '---|' * len(C.ARMS1) + '\n')
    for i, b in enumerate(BARS):
        fh.write('| %.0e | ' % b + ' | '.join(
            '%.0f' % (prof[n]['day'][i] - prof[n]['night'][i]) for n, _ in C.ARMS1)
            + ' |\n')
    fh.write('\n### differences from prodbin > 50 K\n\n')
    for reg in ['day', 'term', 'night']:
        for i, b in enumerate(BARS):
            for n, _ in C.ARMS1[1:]:
                dT = prof[n][reg][i] - prof['prodbin'][reg][i]
                if abs(dT) > 50:
                    fh.write('- %s, %s, %.0e bar: %+.0f K\n' % (n, reg, b, dT))
    fh.write('\n### horizontally averaged T [K] in the eta-capped shell '
             '(r/Rp 1.20-1.28)\n\n| r/Rp | '
             + ' | '.join('%s %s' % (n, r) for n, _ in C.ARMS1 for r in ['night'])
             + ' |\n|---|' + '---|' * len(C.ARMS1) + '\n')
    r = prof['prodbin']['r']
    sel = np.where((r > 1.20) & (r < 1.28))[0][::2]
    for i in sel:
        fh.write('| %.3f | ' % r[i] + ' | '.join(
            '%.0f' % prof[n]['Tr_night'][i] for n, _ in C.ARMS1) + ' |\n')
    fh.close()

    fig, ax = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    col = dict(prodbin='k', sph='C0', sphbeam='C3', noang='C2', offfix='C1')
    for a_, reg, t in zip(ax, ['day', 'term', 'night'],
                          ['day (|lon|<30)', 'terminator', 'night (|lon|>150)']):
        for n, _ in C.ARMS1:
            a_.plot(prof[n][reg], BARS, '-o', ms=3, color=col[n], label=n)
        a_.set_yscale('log')
        a_.invert_yaxis()
        a_.set_xlabel('T [K]')
        a_.set_title(t)
        a_.grid(alpha=.3)
    ax[0].set_ylabel('p [bar]')
    ax[0].legend(fontsize=8)
    fig.suptitle('T(p) at the last dump of each arm (rot ~283.4)')
    fig.tight_layout()
    fig.savefig(C.OUT + 'fig_tp.png', dpi=130)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for n, _ in C.ARMS1:
        ax[0].plot(prof[n]['r'], prof[n]['Tr_night'], color=col[n], label=n + ' night')
        ax[0].plot(prof[n]['r'], prof[n]['Tr_day'], '--', color=col[n], alpha=.6)
        ax[1].plot(prof[n]['r'], prof[n]['Tr_night'] - prof['prodbin']['Tr_night'],
                   color=col[n], label=n)
    for a_ in ax:
        a_.axvspan(1.230, 1.245, color='orange', alpha=.2)
        a_.set_xlabel('r / Rp')
        a_.grid(alpha=.3)
        a_.set_xlim(1.0, 1.6)
    ax[0].set_ylabel('horizontally averaged T [K]  (solid night, dashed day)')
    ax[0].legend(fontsize=8)
    ax[1].set_ylabel('night-side T - prodbin [K]')
    ax[1].set_title('shaded: the shell that sets dt in the sph arms')
    fig.tight_layout()
    fig.savefig(C.OUT + 'fig_Tr.png', dpi=130)
    print(open(C.OUT + 'tp_tables.md').read())


main()
