import sys, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/src'); import eos_gen as eg
A = 7.5657e-15
def tot(rho, T, comp):                 # gas (code EOS) + radiation: specific e, P
    es, ps, mu, xe = eg.evaluate(rho, T, **comp)
    return es + A*T**4/rho, ps*rho + A*T**4/3
def gamma1(rho, T, comp, d=1e-4):
    e0, p0 = tot(rho, T, comp)
    eT = (tot(rho, T*(1+d), comp)[0] - tot(rho, T*(1-d), comp)[0]) / (2*d*T)
    er = (tot(rho*(1+d), T, comp)[0] - tot(rho*(1-d), T, comp)[0]) / (2*d*rho)
    pT = (np.log(tot(rho, T*(1+d), comp)[1]) - np.log(tot(rho, T*(1-d), comp)[1])) / (2*d)
    pr = (np.log(tot(rho*(1+d), T, comp)[1]) - np.log(tot(rho*(1-d), T, comp)[1])) / (2*d)
    dTdrho_s = (p0/rho**2 - er) / eT          # first law at constant s
    return pr + pT * dTdrho_s * rho / T, p0
for name, f, comp, GM, R in [
  ('RSG E', '/resnick/groups/carnegie_poc/jingze/rsgE_page/ic/ic_rsg_g25_deep.txt', dict(xh=0.653848, yhe=0.326395, a_metal=16.0), 1.48422650e27, 5.80019e13),
  ('HG', '/resnick/groups/carnegie_poc/jingze/hgglob_1021/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt', dict(xh=0.0009496, yhe=0.97752, a_metal=1.8047), 3.61091955e26, 7.94653e12)]:
    a = np.loadtxt(f, comments='#'); r, rho, T = a[:, 0], a[:, 1], a[:, 5]
    g1, P = gamma1(rho, T, comp); g = GM / r**2
    N2 = g * (np.gradient(np.log(P), r) / g1 - np.gradient(np.log(rho), r))
    Hp = P / (rho * g); sig = np.sqrt(np.clip(-N2, 0, None)); x = r / R
    m = x < 1.0; i = np.argmax(np.where(m, sig, 0)); smax = sig[i]
    print(f'== {name}: unstable (N2<0) in r/Rph {x[m & (N2 < 0)].min():.3f}-{x[m & (N2 < 0)].max():.4f}')
    print(f'   max growth rate at r/Rph {x[i]:.4f}: 1/sigma {1/smax/86400:.3f} d, H_p/R {Hp[i]/R:.4f}, l ~ pi r/H_p {np.pi*r[i]/Hp[i]:.0f}, Gamma1 {g1[i]:.3f}')
    for fr in [0.5, 0.1, 0.01]:
        w = m & (sig >= fr*smax); print(f'   sigma >= {fr:4.2f} max: r/Rph {x[w].min():.4f}-{x[w].max():.4f}  (r {r[w].min():.4e}-{r[w].max():.4e} cm), l at bottom {np.pi*r[w].min()/Hp[w][0]:.0f}')
    for q in [0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998]:
        j = np.argmin(abs(x - q)); print(f'   r/Rph {q:5.3f}: 1/sigma {1/max(sig[j],1e-30)/86400:10.3f} d  H_p/R {Hp[j]/R:.4f}  l~{np.pi*r[j]/Hp[j]:5.0f}  G1 {g1[j]:.3f}')
