#!/usr/bin/env python3
"""Is the table run's growing supersonic atmosphere caused by the fixed-T top lid?

    python3 lidtest.py [run ...]        default: table lidtest/lid4860 lidtest/lid3000
                                                 lidtest/lid2500

Four 96x64^2 runs identical except for problem/T_top_fix (3500 K = the `table` baseline,
tuned for the IDEAL-gas run; 4860 K = the IC's own temperature minimum; 3000 and 2500 K =
closer to what the tabulated-EOS atmosphere actually settles to). For each it reports the
supersonic fraction and the top-of-box structure versus time, and writes

    soleos/plots/lidtest.png       f(M>1)(t), KE(t), and <T>(z) near the lid
    soleos/lidtest.txt             the numbers

Reading the result:
  * if the growth of f(M>1) tracks T_top_fix, the lid temperature drives it -- the lid is
    inconsistent with the atmosphere the general EOS wants and is forcing the top;
  * if all four grow at much the same rate, the lid TEMPERATURE is not the driver. The
    boundary pins the ghost internal energy at any value of T_top_fix, so it reflects
    the acoustic waves convection sends up, and the fix is a damping layer or a taller
    box rather than a different number.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(_HERE, os.pardir, os.pardir, 'vis', 'python'))
sys.path.insert(0, _HERE)
import compare                                        # noqa: E402
import eoslib                                         # noqa: E402

ROOT = compare.ROOT
DEFAULT = ['table', 'lidtest/lid4860', 'lidtest/lid3000', 'lidtest/lid2500']


def series(run, eos):
    """Supersonic fraction, KE and the near-lid temperature profile at every snapshot."""
    import dumps
    files = dumps.dump_files(os.path.join(ROOT, run))
    out = {'t': [], 'fsup': [], 'machmax': [], 'Ttop': [], 'Tmin': []}
    for i, f in enumerate(files):
        if f is None:
            continue
        d = compare.load(run, i)
        rho = d['dens']
        espec = d['eint']/rho
        T = eos.temperature(rho, espec)
        p = eos.pressure(rho, espec)
        cs = np.sqrt(eos.gamma1(rho, espec)*p/rho)
        M = np.sqrt(d['velx']**2 + d['vely']**2 + d['velz']**2)/cs
        out['t'].append(d['time'])
        out['fsup'].append(float((M > 1).mean()))
        out['machmax'].append(float(np.sqrt((M**2).mean((0, 1))).max()))
        out['Ttop'].append(float(T.mean((0, 1))[-1]))
        out['Tmin'].append(float(T.mean((0, 1)).min()))
        out['Tprof'] = T.mean((0, 1))          # last snapshot only
        out['z'] = d['z']/d['z'][-1]
    for k in ['t', 'fsup', 'machmax', 'Ttop', 'Tmin']:
        out[k] = np.array(out[k])
    return out


def main():
    runs = sys.argv[1:] or DEFAULT
    # tag the outputs by run set, so a second comparison does not overwrite the first
    tag = '' if runs == DEFAULT else '_' + '_'.join(r.split('/')[-1] for r in runs)
    res = {}
    lines = []

    def w(s=''):
        lines.append(s)
        print(s)

    for r in runs:
        res[r] = series(r, eoslib.eos_for_run(os.path.join(ROOT, r), compare.GRID))

    w('lid test: does problem/T_top_fix drive the supersonic atmosphere?')
    w('=' * 88)
    for r in runs:
        s = res[r]
        w('')
        w('%s' % r)
        w('    t      f(M>1)   max Mach_rms   <T> at lid   min <T>')
        for i in range(len(s['t'])):
            w('%7.0f   %6.2f%%      %5.2f         %6.0f K    %6.0f K'
              % (s['t'][i], 100*s['fsup'][i], s['machmax'][i], s['Ttop'][i],
                 s['Tmin'][i]))

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
    for r in runs:
        s = res[r]
        axs[0].plot(s['t'], 100*s['fsup'], marker='o', ms=3, label=r)
        axs[1].plot(s['t'], s['machmax'], marker='o', ms=3, label=r)
        axs[2].plot(s['Tprof'], s['z'], label=r)
    axs[0].set_ylabel(r'cells with Mach $>1$ [%]')
    axs[1].set_ylabel(r'max over $z$ of Mach$_{\rm rms}$')
    for ax in axs[:2]:
        ax.set_xlabel('t [s]')
    axs[1].axhline(1.0, color='k', lw=0.8)
    axs[2].set_xlabel(r'$\langle T\rangle$ [K], last snapshot')
    axs[2].set_ylabel(r'$z\ /\ z_{\rm max}$')
    axs[2].set_ylim(0.6, 1.0)
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle('solar_convection, tabulated EOS: sensitivity to the fixed-T top lid')
    fig.tight_layout()
    out = os.path.join(ROOT, 'plots', 'lidtest%s.png' % tag)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    w('')
    w('plot in %s' % out)
    with open(os.path.join(ROOT, 'lidtest%s.txt' % tag), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
