"""growth times vs dynamical times at the FT density; r_seed summary. usage: venv/bin/python tgrow.py"""
import pickle, numpy as np, grainlib as gl
ch = gl.Chem()
print('| star | r/R | rho_FT | n_H | t_0.1um a=1 Tg=400/1700 [s] | t_0.3um a=1 Tg=400/1700 [s] | t_ff=sqrt(r^3/GM) [s] |')
print('|---|---|---|---|---|---|---|')
for sn in gl.STARS:
    st = gl.Star(sn)
    for x in [1.5, 2, 3, 4, 5]:
        r = x * st.R; rho = st.rho_ft(r); nH = rho / ch.mH
        t = lambda a, Tg: a / (gl.V_MON * nH * ch.eMg * np.sqrt(8 * gl.kB * Tg / (np.pi * gl.M_MG)) / 8)
        print('| %s | %.1f | %.1e | %.1e | %.1e / %.1e | %.1e / %.1e | %.1e |' % (
            sn, x, rho, nH, t(1e-5, 400), t(1e-5, 1700), t(3e-5, 400), t(3e-5, 1700),
            np.sqrt(r**3 / (gl.G * st.M))))
res = pickle.load(open(gl.D + 'grid.pkl', 'rb'))
print('\nr_seed/R (first growth), dens=ft, median [min-max] over runs that formed grains:')
for sn in gl.STARS:
    for var in ['lowk', 'pure', 'Fe3e-4', 'Fe1e-3']:
        rs = np.array([r['r_seed'] for r in res if r['star'] == sn and r['var'] == var
                       and r['dens'] == 'ft' and np.isfinite(r['r_seed'])])
        print('  %s %-6s n=%5d  %.2f [%.2f-%.2f]' % (sn, var, len(rs), np.median(rs), rs.min(), rs.max()))
# drift and vinf over dust-driven winds
dw = [r for r in res if r['outcome'] == 'wind' and r['rmax0'] < 20 and r['vinf'] > 0]
for sn in gl.STARS:
    d = [r for r in dw if r['star'] == sn]
    v = np.array([r['vinf'] for r in d]) / 1e5; vd = np.array([r['vdmax'] for r in d]) / 1e5
    md = np.array([r['mdot'] for r in d]) / gl.Msun * gl.YR
    print('%s dust-driven winds n=%d: vinf median %.0f [%.0f-%.0f] km/s; vd_max median %.0f [%.0f-%.0f]; '
          'Mdot(>v0) range %.1e-%.1e' % (sn, len(d), np.median(v), v.min(), v.max(), np.median(vd),
                                          vd.min(), vd.max(), md.min(), md.max()))
for sn in gl.STARS:
    st = gl.Star(sn)
    for gg in [0, 0.1, 0.3]:
        ve = st.vesc * np.sqrt(1 - gg)
        print('%s Ggas=%.1f: ballistic escape v0 > %.1f km/s -> Mdot(>v_esc,eff) = %.1e Msun/yr' %
              (sn, gg, ve / 1e5, st.mdot_gt(ve) / gl.Msun * gl.YR))
