"""arms.py NAME:NPZ:HST ... : RSGARMS read-out (conv_on.py npz + user.hst): KE_lat, KE_r (780-836 Rsun), v_lat rms at 815/825 Rsun, L_top/L_in (5 d mean); table + figure"""
import sys, re, numpy as np
DAY, RS = 86400., 6.957e10
def hst(p):
    h = [l for l in open(p) if l.startswith('#') and '[' in l][-1]; n = [m.group(2) for m in re.finditer(r'\[(\d+)\]=(\S+)', h)]
    a = np.loadtxt(p); a = a[np.unique(a[:, 0], return_index=True)[1]]; return a[:, 0]/DAY, a[:, n.index('L_top')]/a[:, n.index('L_in')]
R = {}
for arg in sys.argv[1:]:
    nm, f, hp = arg.split(':'); S = np.load(f); r = S['r']; dV = 4*np.pi*r**2*np.diff(S['xf']); lay = (r > 780*RS) & (r < 836*RS)
    t = S['t']/DAY; th, lt = hst(hp)
    d = dict(t=t, kel=(S['ke_l']*dV*lay).sum(1), ker=(S['ke_r']*dV*lay).sum(1), v815=S['vl_rms'][:, np.abs(r - 815*RS).argmin()]/1e5,
             v825=S['vl_rms'][:, np.abs(r - 825*RS).argmin()]/1e5, lt=np.array([lt[(th > x - 5) & (th <= x)].mean() if ((th > x - 5) & (th <= x)).any() else np.nan for x in t]), th=th, lth=lt)
    R[nm] = d
print('%-8s %6s %10s %10s %7s %7s %7s' % ('run', 't[d]', 'KE_lat', 'KE_r', 'vl815', 'vl825', 'Ltop/Lin'))
for nm, d in R.items():
    for k in range(len(d['t'])):
        if d['t'][k] >= 29: print('%-8s %6.1f %10.3e %10.3e %7.2f %7.2f %7.4f' % (nm, d['t'][k], d['kel'][k], d['ker'][k], d['v815'][k], d['v825'][k], d['lt'][k]))
if len(sys.argv) > 1 and __name__ == '__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    import os; sys.path.insert(0, os.environ.get('PLOTSTYLE', '/orion/ptmp/jinma/clean_1115/stag2/plot_style')); import idort_style; idort_style.use()
    fig, ax = plt.subplots(4, 1, figsize=(10, 14), sharex=True); fig.subplots_adjust(hspace=0.01)
    cl = dict(oldS='k', ctl='0.5', xmono='r', tanoff='b', newS='C2', nrigid='C2', oldSam='C1', oldS2='k', OC='k', OA='r', NC='C2', NI='b')
    for nm, d in R.items():
        c = cl.get(nm, 'C4'); ls = '--' if nm in ('newS', 'oldS2', 'NC', 'NI') else ('-.' if nm == 'nrigid' else '-')
        ax[0].semilogy(d['t'], d['kel'], ls, color=c, marker='x', label=nm); ax[1].semilogy(d['t'], d['ker'], ls, color=c, marker='x')
        ax[2].plot(d['t'], d['v815'], ls, color=c, marker='x'); ax[2].plot(d['t'], d['v825'], ':', color=c)
        k = d['th'] > 1; ax[3].plot(d['th'][k][::20], d['lth'][k][::20], ls, color=c, lw=1)
    ax[0].set_ylabel(r'$\rm KE_{lat}\ [erg]$'); ax[1].set_ylabel(r'$\rm KE_{r}\ [erg]$'); ax[2].set_ylabel(r'$\rm v_{lat}\ [km\,s^{-1}]$'); ax[3].set_ylabel(r'$\rm L_{top}/L_{in}$')
    ax[2].text(0.02, 0.85, r'$\rm 815\,R_\odot$ (solid), $\rm 825\,R_\odot$ (dotted)', transform=ax[2].transAxes, fontsize=16)
    ax[3].set_xlabel(r'$\rm t\ [d]$'); ax[3].set_xlim(25, 110); ax[3].set_ylim(0.35, 0.95); ax[0].legend(fontsize=18, frameon=False, ncol=2)
    fig.savefig(os.environ.get('OUT', 'arms.png'), dpi=300, bbox_inches='tight'); print('arms.png')
