"""KE_int(t) of the page's run ($BSG_REPRO_RUN, default repro_sf_4n: 1 % T seed,
force_reference_work = split) against the unseeded force_reference_work = full run
($BSG_REPRO_FULL_RUN, default repro_4n, stopped), from their user.hst (column 11 KE_int;
read only), with the paper's kinetic-energy level (~1e43 erg, Ma+2026 Fig. 2) as a
reference line.  Writes out/ours_ke.png and cache/ke.json.

usage: python3 ke_compare.py"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bsg_figs as B  # noqa: E402
plt = B.plt

S = os.path.dirname(os.path.abspath(__file__))
PAPER_KE = 1e43


def load(run):
    h = B.read_hst(os.path.join(run, 'bsg3d.user.hst'))
    t = h['time'] / B.DAY
    ke = h['KE_int']
    o = np.argsort(t, kind='stable')
    return t[o], ke[o]


def at(t, y, td):
    i = min(np.searchsorted(t, td), len(t) - 1)
    return float(y[i])


def main():
    rs, rf = B.ARMS['repro']['run'], B.ARMS['repro_full']['run']
    ts, ks = load(rs)
    tf, kf = load(rf)
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    m = kf > 0
    ax.semilogy(tf[m], kf[m], color='0.45', lw=1.4,
                label='unseeded, force work = full (%s, stopped)' % os.path.basename(rf))
    m = ks > 0
    ax.semilogy(ts[m], ks[m], color='C0', lw=1.6,
                label='1 %% T seed, force work = split (%s)' % (
                    'repro_sf_4n, Raven + viper' if 'merged' in rs else os.path.basename(rs)))
    ax.axhline(PAPER_KE, color='C3', ls='--', lw=1.3,
               label=r'Ma+2026 Fig. 2 level ($\approx 10^{43}$ erg)')
    ax.set_xlim(0, max(ts.max(), tf.max()) * 1.02)
    ax.set_ylim(1e38, 3e43)
    ax.set_xlabel('Time since Simulation Starts [days]')
    ax.set_ylabel(r'KE$_{\rm int}$ [erg]')
    ax.grid(True, which='major', alpha=0.3)
    ax.legend(loc='lower right', fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(S, 'out', 'ours_ke.png'), dpi=150)
    plt.close(fig)
    late = ts >= 5.0
    res = dict(run=rs, run_full=rf, t_end_d=float(ts.max()), t_end_full_d=float(tf.max()),
               ke_now=float(ks[-1]), ke_full_7d=at(tf, kf, 7.0), ke_full_end=float(kf[-1]),
               ke_full_max=float(kf.max()), t_full_max_d=float(tf[np.argmax(kf)]),
               ke_seed_max=float(ks.max()), t_seed_max_d=float(ts[np.argmax(ks)]),
               ke_seed_min_after5=float(ks[late].min()) if late.any() else None,
               ke_seed_max_after5=float(ks[late].max()) if late.any() else None,
               ke_seed_7d=at(ts, ks, 7.0), paper_ke=PAPER_KE)
    json.dump(res, open(os.path.join(S, 'cache', 'ke.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
