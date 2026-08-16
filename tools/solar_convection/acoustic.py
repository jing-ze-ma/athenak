#!/usr/bin/env python3
"""Does the atmosphere trap acoustic waves, or just fail to damp them?

    python3 acoustic.py [run ...]          default: ideal table

The tabulated-EOS runs develop a supersonic atmosphere that keeps growing, and the lid
test showed the growth is insensitive to problem/T_top_fix (see soleos/lidtest.txt). Two
explanations survive:

  (1) REFLECTION -- the top boundary pins the ghost internal energy at any T_top_fix, so
      upgoing waves cannot leave and their energy accumulates;
  (2) WEAK DAMPING -- the tabulated atmosphere is ~100x more rarefied, so kappa*rho and
      hence the radiative damping of acoustic waves is far weaker than in the ideal run,
      and the waves steepen into shocks.

Both runs share the SAME boundary condition, so (1) alone cannot explain why only the
table run runs away -- the ideal run is the control that separates them. This script
measures both effects directly.

Method. Split each horizontal plane into a mean and a fluctuation, rho0(z) = <rho>, and
likewise p0, cs0, vz0. For a linear acoustic wave the up- and down-going characteristics
are w' +- p'/(rho0 cs0), carrying the fluxes

    F+ = (rho0 cs0/4) < (w' + p'/(rho0 cs0))^2 >     upward
    F- = (rho0 cs0/4) < (w' - p'/(rho0 cs0))^2 >     downward
    F_net = <p' w'> = F+ - F-                        (exact identity, a useful check)

so R = F-/F+ is a reflection coefficient: R -> 1 is a closed lid (standing wave, no net
escape), R -> 0 is free transmission. The wave energy density is

    E_w = (1/2) rho0 <w'^2> + (1/2) <p'^2>/(rho0 cs0^2)      kinetic + compressional

and the damping is quantified by the optically thin (Spiegel) radiative relaxation time
against the acoustic time on a pressure scale height,

    t_rad = c_v/(16 kappa sigma T^3),    t_ac = H_p/cs0,    H_p = p0/(rho0 g)

with c_v = de/dT taken from the run's OWN eos by finite difference, and kappa the pgen's
own get_kapr. t_rad/t_ac >> 1 means waves cross a scale height untouched by radiation.

Reading the result. Compare the two runs at the same height above the photosphere:
  * if R_top is much closer to 1 in the table run, the boundary is trapping and the fix is
    a sponge layer or a taller box;
  * if R_top is similar in both but t_rad/t_ac is far larger in the table run, reflection
    is a shared condition and the DIFFERENCE is damping -- the growth is then physical,
    and the honest write-up is that the run needs radiative damping (or an explicit
    sponge) to reach a steady atmosphere, not that it is broken.
The energy budget printed at the end is the third leg: if dE_atm/dt is a large fraction of
the acoustic luminosity crossing the photosphere, essentially nothing is escaping or being
dissipated.

Caveat: the decomposition is linear. It is quantitative while the atmosphere is subsonic
(t < ~12000 s here) and only indicative once Mach > 1, so the EARLY times carry the
argument -- which is the point, since that is when the growth starts.

RESULT (2026-08-15, the 96x64^2 pair): NEITHER hypothesis survives. R is ~0.3-0.7 and
indistinguishable between the two runs -- the boundary reflects partially in both -- and
t_rad/t_ac is ~11-17 at mid-atmosphere in BOTH, so both atmospheres are radiatively
undamped. The acoustic luminosity crossing the photosphere is also the same, ~5-10e9 in
both, and steady. What differs is the AMPLITUDE that luminosity reaches: the tabulated
atmosphere is ~10x thinner at mid-atmosphere (rho0 1.1e-8 vs 9.7e-8) with a 34% smaller
c_s, so at t=10000 it carries a SMALLER flux (2.7e7 vs 7.4e7) at a 1.8x LARGER w'_rms and
2.7x larger Mach (0.109 vs 0.040). That is F ~ rho0 cs0 <w'^2> read backwards. The ideal
run stays linear forever; the table run reaches shock amplitude, and shock heating then
inflates the atmosphere, which lets still more flux through -- by t=20000 F+ at
mid-atmosphere has grown 10x, to 3x the ideal run's, at Mach 0.34. The rarefaction itself
is real EOS physics (H2 formation raises mu, halving the scale height). So the runaway is
physical amplification, not a boundary artifact -- but this box is too short and too
undamped to host it, so the atmospheric statistics of these runs must not be quoted as
converged. A taller box with a sponge layer, or real radiative damping, is what it needs.
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
import dumps                                          # noqa: E402
import eoslib                                         # noqa: E402

ROOT = compare.ROOT
GRAV = 2.74e4              # |g|, cm/s^2, as in the pgen
SIGMA = 5.6704e-5          # Stefan-Boltzmann, cgs
DEFAULT = ['ideal', 'table']


def eos_for(run):
    return eoslib.eos_for_run(os.path.join(ROOT, run), compare.GRID)


def cv_of(eos, rho, espec, frac=1.0e-3):
    """Specific heat at constant density, de/dT, by finite difference through the EOS."""
    de = frac*espec
    t0 = eos.temperature(rho, espec)
    t1 = eos.temperature(rho, espec + de)
    return de/np.maximum(t1 - t0, 1.0e-30)


def snapshot(run, eos, idx):
    """Acoustic decomposition of one snapshot, as horizontally averaged z-profiles."""
    d = compare.load(run, idx)
    rho = d['dens']
    espec = d['eint']/rho
    T = eos.temperature(rho, espec)
    p = eos.pressure(rho, espec)
    cs = np.sqrt(eos.gamma1(rho, espec)*p/rho)
    vz = d['velx']                                    # x1 is the vertical axis

    ax = (0, 1)
    rho0, p0, cs0, vz0 = rho.mean(ax), p.mean(ax), cs.mean(ax), vz.mean(ax)
    pp = p - p0                                       # broadcasting over (nx3, nx2, nx1)
    wp = vz - vz0
    q = pp/(rho0*cs0)                                 # p'/(rho0 cs0), a velocity

    wrms = np.sqrt((wp**2).mean(ax))
    fplus = 0.25*rho0*cs0*((wp + q)**2).mean(ax)
    fminus = 0.25*rho0*cs0*((wp - q)**2).mean(ax)
    fnet = (pp*wp).mean(ax)
    ew = 0.5*rho0*(wp**2).mean(ax) + 0.5*(pp**2).mean(ax)/(rho0*cs0**2)

    # radiative relaxation vs acoustic crossing time
    cv = cv_of(eos, rho, espec).mean(ax)
    kap = eoslib.get_kapr(rho, T).mean(ax)
    T0 = T.mean(ax)
    t_rad = cv/(16.0*kap*SIGMA*T0**3)
    hp = p0/(rho0*GRAV)
    t_ac = hp/cs0

    # optical depth from the top, to locate the photosphere
    tau = np.cumsum((kap*rho0*d['dz'])[::-1])[::-1]
    return {'t': d['time'], 'z': d['z'], 'dz': d['dz'], 'rho0': rho0, 'cs0': cs0,
            'fplus': fplus, 'fminus': fminus, 'fnet': fnet, 'ew': ew, 'wrms': wrms,
            'R': fminus/np.maximum(fplus, 1.0e-300), 't_rad': t_rad, 't_ac': t_ac,
            'tau': tau, 'T0': T0}


def series(run):
    eos = eos_for(run)
    files = dumps.dump_files(os.path.join(ROOT, run))
    out = [snapshot(run, eos, i) for i, f in enumerate(files) if f is not None]
    # the photosphere of the LAST snapshot fixes the base of "the atmosphere"
    tau = out[-1]['tau']
    kph = int(np.argmin(np.abs(tau - 2.0/3.0)))
    for s in out:
        s['kph'] = kph
        # wave energy stored above the photosphere, per unit area
        s['E_atm'] = float(np.sum(s['ew'][kph:])*s['dz'])
        s['F_ph'] = float(s['fplus'][kph])            # acoustic luminosity going up
        s['R_top'] = float(s['R'][-1])
        s['R_mid'] = float(s['R'][(kph + len(tau))//2])
    return out


def main():
    runs = sys.argv[1:] or DEFAULT
    # tag the outputs by run set, so a second comparison does not overwrite the first
    tag = '' if runs == DEFAULT else '_' + '_'.join(r.split('/')[-1] for r in runs)
    res = {r: series(r) for r in runs}
    lines = []

    def w(s=''):
        lines.append(s)
        print(s)

    w('acoustic budget of the atmosphere: reflection vs damping')
    w('=' * 92)
    for r in runs:
        s = res[r]
        z = s[-1]['z']/s[-1]['z'][-1]
        kph = s[-1]['kph']
        w('')
        w('%s   (photosphere at z/zmax = %.3f, index %d)' % (r, z[kph], kph))
        w('    t      F+ at ph     E_atm      dE/dt      dE/dt / F+     R mid   R top   '
          'min t_rad/t_ac')
        for i, ss in enumerate(s):
            if i == 0:
                dedt, frac = np.nan, np.nan
            else:
                dt = ss['t'] - s[i-1]['t']
                dedt = (ss['E_atm'] - s[i-1]['E_atm'])/dt
                frac = dedt/ss['F_ph'] if ss['F_ph'] > 0 else np.nan
            ratio = (ss['t_rad']/ss['t_ac'])[kph:]
            w('%7.0f  %10.3e  %10.3e  %+10.3e  %+10.3f     %6.3f  %6.3f   %10.1f'
              % (ss['t'], ss['F_ph'], ss['E_atm'], dedt, frac,
                 ss['R_mid'], ss['R_top'], ratio.min()))

    # ---- figure ------------------------------------------------------------------------
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for n, r in enumerate(runs):
        s = res[r]
        c = colors[n % len(colors)]
        last, early = s[-1], s[min(10, len(s)-1)]
        zz = last['z']/last['z'][-1]
        kph = last['kph']
        axs[0, 0].plot(last['fplus'], zz, c, ls='-', label='%s  $F_+$' % r)
        axs[0, 0].plot(last['fminus'], zz, c, ls='--', label='%s  $F_-$' % r)
        axs[0, 1].plot(early['R'], zz, c, ls='--', label='%s  t=%.0f' % (r, early['t']))
        axs[0, 1].plot(last['R'], zz, c, ls='-', label='%s  t=%.0f' % (r, last['t']))
        axs[0, 2].plot(last['t_rad']/last['t_ac'], zz, c, label=r)
        axs[1, 0].plot([ss['t'] for ss in s], [ss['E_atm'] for ss in s], c,
                       marker='o', ms=3, label=r)
        axs[1, 1].plot([ss['t'] for ss in s], [ss['R_top'] for ss in s], c,
                       marker='o', ms=3, label='%s  top' % r)
        axs[1, 1].plot([ss['t'] for ss in s], [ss['R_mid'] for ss in s], c, ls='--',
                       marker='s', ms=3, label='%s  mid-atm' % r)
        axs[1, 2].plot([ss['t'] for ss in s], [ss['F_ph'] for ss in s], c,
                       marker='o', ms=3, label=r)
        for ax in [axs[0, 0], axs[0, 1], axs[0, 2]]:
            ax.axhline(zz[kph], color=c, lw=0.7, ls=':')

    axs[0, 0].set_xscale('log')
    axs[0, 0].set_xlabel(r'$F_\pm$ [erg cm$^{-2}$ s$^{-1}$], last snapshot')
    axs[0, 1].set_xlabel(r'$R = F_-/F_+$   (1 = closed lid, 0 = free escape)')
    axs[0, 1].set_xlim(0, 1.4)
    axs[0, 1].axvline(1.0, color='k', lw=0.8)
    axs[0, 2].set_xscale('log')
    axs[0, 2].set_xlabel(r'$t_{\rm rad}/t_{\rm ac}$  ($\gg 1$: no radiative damping)')
    axs[0, 2].axvline(1.0, color='k', lw=0.8)
    for ax in axs[0]:
        ax.set_ylabel(r'$z\ /\ z_{\rm max}$')
    axs[1, 0].set_ylabel(r'wave energy above the photosphere [erg cm$^{-2}$]')
    axs[1, 1].set_ylabel(r'$R = F_-/F_+$')
    axs[1, 1].axhline(1.0, color='k', lw=0.8)
    axs[1, 2].set_ylabel(r'$F_+$ at the photosphere [erg cm$^{-2}$ s$^{-1}$]')
    axs[1, 0].set_yscale('log')
    for ax in axs[1]:
        ax.set_xlabel('t [s]')
    for ax in axs.flat:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle('solar_convection: is the atmosphere trapping waves, or failing to '
                 'damp them?')
    fig.tight_layout()
    out = os.path.join(ROOT, 'plots', 'acoustic%s.png' % tag)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    # ---- the amplitude argument, at a common height ------------------------------------
    w('')
    w('=' * 92)
    w('mid-atmosphere state (halfway from the photosphere to the lid). The SAME acoustic')
    w('luminosity driven into a thinner atmosphere means a larger velocity amplitude:')
    w('F+ ~ rho0 cs0 <w\'^2>, so w\'_rms ~ sqrt(F+/(rho0 cs0)).')
    w('')
    w('%-22s %10s %10s %10s %10s %10s %10s' % ('run / t', 'rho0', 'cs0 [km/s]',
                                               "w'rms", 'Mach_rms', 'F+', 't_rad/t_ac'))
    w('-' * 92)
    for r in runs:
        for ss in [res[r][min(10, len(res[r])-1)], res[r][-1]]:
            k = (ss['kph'] + len(ss['z']))//2
            w('%-22s %10.3e %10.2f %10.3e %10.3f %10.3e %10.2f'
              % ('%s  t=%.0f' % (r, ss['t']), ss['rho0'][k], ss['cs0'][k]/1e5,
                 ss['wrms'][k], ss['wrms'][k]/ss['cs0'][k], ss['fplus'][k],
                 (ss['t_rad']/ss['t_ac'])[k]))
    w()
    w('plot in %s' % out)
    with open(os.path.join(ROOT, 'acoustic%s.txt' % tag), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
