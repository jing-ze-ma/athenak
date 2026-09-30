import sys, numpy as np
sys.path.insert(0, 'wt/docs/handover/scripts')
import dhj_remap as d
h = d.read_rst('/viper/ptmp2/jinma/w121prod_0929/w1x/rst/dhj.00600.rst')
par = h['par']; ng = h['ng']
nx1o, r0, r1, co = d.grid_of(par)
eo = d.edges(nx1o, r0, r1, co, ng); ro = d.centroids(eo)
g = float(par['problem']['grav']); apl = float(par['problem']['ap'])
om2 = float(par['problem']['omega'])**2
phig = g*apl*(1-apl/ro)
full = dict(rho=h['u'][:, 0], m1=h['u'][:, 1], m2=h['u'][:, 2], m3=h['u'][:, 3],
            ei=h['eint'], p=h['p'], E=h['u'][:, 4])
PA = {k: d.panel_arrays(h, full[k], 32, 32) for k in full}
PU = {k: PA[k]/PA['rho'] for k in ('m1', 'm2', 'm3')}
f = 8
rows = []
for pn in range(6):
    hp = d.HProlong(f, 0, 32, 0, 32, 32, 32, pn)
    ch = {}
    for k in ('rho', 'ei'):
        c, sj, sk = d.slopes(PA[k], pn, 0, 32, 0, 32)
        ch[k], _ = hp.children(c, sj, sk, positive=True, qmin=1e-16 if k == 'rho' else 0)
    rW = hp.parent_sum(ch['rho'])
    W = hp.w4.sum(axis=(1, 3))[..., None]
    for k in ('m1', 'm2', 'm3'):
        c, sj, sk = d.slopes(PU[k], pn, 0, 32, 0, 32)
        uc, _ = hp.children(c, sj, sk)
        mp = PA[k][pn, 1:33, 1:33]*W
        du = (mp - hp.parent_sum(ch['rho']*uc))/rW
        ch[k] = ch['rho']*(uc.reshape(32, f, 32, f, -1) + du[:, None, :, None, :]).reshape(256, 256, -1)
    kef = d.kinetic(ch['rho'], ch['m1'], ch['m2'], ch['m3'], hp.cosc[..., None])
    B_f = -0.5*om2*hp.sin2
    rpf = ch['rho']*(phig + B_f[..., None]*ro**2)
    pc = {k: PA[k][pn, 1:33, 1:33] for k in PA}
    # parent-level quantities
    w, cosp, sin2p = d.panel_cols(0, 32, 0, 32, 32, 32, pn)
    kep = d.kinetic(pc['rho'], pc['m1'], pc['m2'], pc['m3'], cosp[..., None])
    rpp = pc['rho']*(phig + (-0.5*om2*sin2p)[..., None]*ro**2)
    eW = pc['ei']*W
    dK = (hp.parent_sum(kef) - kep*W)/eW
    dP = (hp.parent_sum(rpf) - rpp*W)/eW
    dC = (pc['E'] - pc['ei'] - kep - rpp)/pc['ei']
    sel = np.zeros_like(dK, bool); sel[..., ng:ng+nx1o] = True
    ff = 1 - dK - dP + dC
    bad = sel & (np.abs(ff-1) > 0.1)
    for idx in zip(*np.nonzero(bad)):
        rows.append((pc['p'][idx]/1e6, dK[idx], dP[idx], dC[idx], idx[2]-ng, pn))
rows = np.array(rows)
print(len(rows))
o = np.argsort(-rows[:, 0])
np.set_printoptions(precision=3, linewidth=150)
print('p_bar dK/e dPhi/e dC2P/e i panel'); print(rows[o[:15]])
for lo in (1e-9,1e-6,1e-4,1e-2):
    s = rows[:,0] >= lo; print(lo, s.sum(), 'median|dK| %.3f median|dP| %.3f |dC| %.3f' % (np.median(np.abs(rows[s,1])) if s.any() else 0, np.median(np.abs(rows[s,2])) if s.any() else 0, np.median(np.abs(rows[s,3])) if s.any() else 0))
