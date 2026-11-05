import numpy as np, glob, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
EU = 1e-9*(6.957e10)**3*(6.957e10/1.899769e4)**2
SR = 4*np.pi/((np.cos(np.pi/3)-np.cos(2*np.pi/3))*(np.pi/3))
S = EU*SR/1e45
runs = [('L_Bd1','Bd1 closed flux wall (af 1)','#1f77b4'), ('L_tc4','his BC, taucell 4 (af 0)','#d62728'),
        ('C_Bd1','Bd1 with converged iterations (err 1e-11)','#2ca02c')]
fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
for d, lab, c in runs:
    try: a = np.loadtxt('runs/%s/out/ledger.txt' % d)
    except Exception: continue
    a = a[np.all(np.isfinite(a), axis=1)]; t = a[:, 0]/4.547922
    b = 1 + 2*12   # whole domain
    dE = (a[:, b]-a[0, b])*S
    T = (a[:, b+1:b+12]-a[0, b+1:b+12])*S
    ax[0].plot(t, dE, color=c, lw=2, label=lab + ': measured dE_tot')
    ax[0].plot(t, T[:, 3]-T[:, 4], color=c, lw=1, ls='--', label='   non-telescoping rad transport')
    ax[0].plot(t, T[:, 4], color=c, lw=1, ls=':', label='   consistent rad transport (= boundary rad flux)')
    ax[0].plot(t, T[:, 1]+T[:, 2]+T[:, 0], color=c, lw=0.8, ls='-.', label='   hydro + gravity')
ax[0].set_xlabel('t [d]'); ax[0].set_ylabel('cumulative energy change [1e45 erg, full sphere]'); ax[0].set_title('whole domain: energy by term')
ax[0].legend(fontsize=6.5); ax[0].grid(alpha=0.3); ax[0].set_xlim(0, 30)
for k, (d, lab, c) in enumerate(runs[:2]):
    fs = sorted(glob.glob('runs/%s/out/ledger_prof_*.txt' % d)); f0 = fs[0]
    f5 = [f for f in fs if '022.7' in f][0]
    p0 = np.loadtxt(f0); p5 = np.loadtxt(f5)
    r = p5[:, 1]; dE = np.cumsum(p5[:, 2]-p0[:, 2])*S
    T = (p5[:, 3:]-p0[:, 3:])*S
    a = ax[1+k]
    a.plot(r, dE, 'k', lw=2, label='measured dE_tot(<r)')
    a.plot(r, np.cumsum(T[:, 3]-T[:, 4]), color='#1f77b4', label='non-telescoping rad transport')
    a.plot(r, np.cumsum(T[:, 4]), color='#9467bd', label='consistent rad transport')
    a.plot(r, np.cumsum(T[:, 0]), color='#2ca02c', label='hydro flux')
    a.plot(r, np.cumsum(T[:, 1]+T[:, 2]), color='#8c564b', label='gravity (mass flux x Phi + point-mass source)')
    a.plot(r, np.cumsum(T[:, 6]+T[:, 7]+T[:, 8]), color='#7f7f7f', label='resets + floors')
    a.set_xlim(20, 35); a.set_xlabel('r [Rsun]'); a.set_ylabel('cumulative over r < r [1e45 erg]')
    a.set_title('%s: by term at 5 d, integrated inward-out' % lab); a.grid(alpha=0.3); a.legend(fontsize=7)
fig.suptitle('Athena++ energy ledger (1-D Linhao BSG): where the energy goes')
fig.tight_layout(); fig.savefig('ledger_user.png', dpi=120)
