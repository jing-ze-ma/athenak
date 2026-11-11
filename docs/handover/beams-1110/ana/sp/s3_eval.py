"""Stage 3 exact-solution evaluators (NOTE-sptests 0.4).
usage: s3_eval.py <mode> <case> <arm> [levels...]
  thin  T6a: black-body sphere r_in = 1, flux 1, thin medium: E_ex = 2 (1 - sqrt(1 - 1/r^2)),
        F_ex = 1/r^2 (c = 1).  Reports max / L1 |E/E_ex - 1|, E/E_ex at r = 1.05, 2, 5, 9.5,
        max |F r^2 - 1|.
  dif   T4: spherical diffusion of E0 = e_out + exp(-40 (r - r0)^2) (u = r E obeys the 1-D
        heat equation, D = c / (3 rho kappa_s); the pulse never reaches the ends).  L1 of
        (E - E_ex) over L1 of (E_ex - e_out) at each output time.
  relax T1: 0-D thermal relaxation, interior cells r < 1.8: max rel error of E and T at each
        output time against the ODE (scipy, rtol 1e-12).
Orders between successive levels are printed for every metric."""
import math
import os
import sys

import numpy as np

import spl
import spt


def rcen(case, n):
    c = spt.case_cfg(case)
    if c.get('tref'):
        # a T ladder: the level is the dt divisor on the fixed grid (spt.mk)
        n = 64 if c['test'] == 'relax' else 128
    rf = spl.rfaces(c['r'][0], c['r'][1], n, c['poly'])
    return spl.rcent(rf), c


def run_dir(case, arm, n):
    return os.path.join(spl.RUNS, case, '%s_%d' % (arm, n))


def ev_thin(case, arm, n):
    d = run_dir(case, arm, n)
    t = spl.tab(spl.tabs(d)[-1])
    r, _ = rcen(case, n)
    ex = 2.0*(1.0 - np.sqrt(1.0 - 1.0/r**2))
    rel = t['m1_e']/ex - 1.0
    fr = t['m1_f1']*r*r - 1.0
    pick = ' '.join('%.2f:%+.3f' % (r[i], rel[i]) for i in
                    [np.argmin(np.abs(r - x)) for x in (1.05, 2.0, 5.0, 9.5)])
    st = spl.logstats(d)
    return {'L1': np.mean(np.abs(rel)), 'max': np.max(np.abs(rel)), 'F': np.max(np.abs(fr)),
            'txt': 'cyc %d Picard %.2f | E/E_ex-1 L1 %.2e max %.2e | at r %s | max|Fr2-1| %.2e'
            % (t['cycle'], st.get('pic', -1), np.mean(np.abs(rel)), np.max(np.abs(rel)), pick,
               np.max(np.abs(fr)))}


def dif_exact(r, t, r0, eo, D):
    s = np.linspace(r0 - 1.5, r0 + 1.5, 6001)
    g = s*np.exp(-40.0*(s - r0)**2)
    w = np.gradient(s)
    k = np.exp(-(r[:, None] - s[None, :])**2/(4.0*D*t))/math.sqrt(4.0*math.pi*D*t)
    return eo + (k @ (g*w))/r


def ev_dif(case, arm, n):
    d = run_dir(case, arm, n)
    r, c = rcen(case, n)
    r0 = float(c['prob']['r0'])
    eo = float(c['prob']['e_out'])
    D = float(c['rad']['c_light'])/(3.0*float(c['rad']['kappa_s']))
    out = []
    for f in spl.tabs(d)[1:]:
        t = spl.tab(f)
        ex = dif_exact(r, t['time'], r0, eo, D)
        out.append((math.sqrt(D*t['time']), np.sum(np.abs(t['m1_e'] - ex))/np.sum(ex - eo)))
    st = spl.logstats(d)
    return {'L1': out[-1][1], 'txt': 'Picard %.2f | ' % st.get('pic', -1) + ' '.join(
        'sqrt(Dt) %.3f L1 %.2e' % o for o in out)}


def relax_ode(T0, E0, rho, gm1, a, c, kap, times):
    from scipy.integrate import solve_ivp

    def f(_, y):
        q = c*kap*(a*y[1]**4 - y[0])
        return [q, -q*gm1/rho]
    s = solve_ivp(f, (0.0, max(times)), [E0, T0], t_eval=sorted(times), rtol=1e-12,
                  atol=1e-14, method='Radau')
    return {tt: (s.y[0][i], s.y[1][i]) for i, tt in enumerate(sorted(times))}


def ev_relax(case, arm, n):
    d = run_dir(case, arm, n)
    r, c = rcen(case, n)
    T0, E0 = float(c['prob']['gas_temp']), float(c['prob']['e_out'])
    rho = float(c['prob']['gas_rho'])
    cl, kap = float(c['rad']['c_light']), float(c['rad']['kappa_p'])
    gm1 = 2.0/3.0
    fm, fh = spl.tabs(d)[1:], spl.tabs(d, 'hydro_w')[1:]
    times = [spl.tab(f)['time'] for f in fm]
    ref = relax_ode(T0, E0, rho, gm1, 1.0, cl, kap, times)
    m = r < 1.8
    out = []
    for a, b in zip(fm, fh):
        tm, th = spl.tab(a), spl.tab(b)
        Ee, Te = ref[tm['time']]
        T = th['eint']*gm1/th['dens']
        out.append((tm['time'], np.max(np.abs(tm['m1_e'][m]/Ee - 1)),
                    np.max(np.abs(T[m]/Te - 1))))
    st = spl.logstats(d)
    return {'L1': out[-1][1], 'T': out[-1][2], 'txt': 'Picard %.2f | ' % st.get('pic', -1) +
            ' '.join('t %.2e E %.2e T %.2e' % o for o in out)}


def core_exact(r, kap=10.0, nmu=4000):
    """IDORT 4.1.2: J, H of a uniform sphere R = 1, chi = kap, B = 1, vacuum outside"""
    x, w = np.polynomial.legendre.leggauss(nmu)
    J = np.zeros_like(r)
    H = np.zeros_like(r)
    for q, rr in enumerate(r):
        d = 1.0 - rr*rr*(1.0 - x*x)
        sq = np.sqrt(np.maximum(d, 0.0))
        s0 = np.maximum(rr*x - sq, 0.0)
        s1 = np.maximum(rr*x + sq, 0.0)
        ell = np.where(d > 0.0, s1 - s0, 0.0)
        inten = 1.0 - np.exp(-kap*ell)
        J[q] = 0.5*np.sum(w*inten)
        H[q] = 0.5*np.sum(w*inten*x)
    return J, H


def ev_core(case, arm, n):
    d = run_dir(case, arm, n)
    t = spl.tab(spl.tabs(d)[-1])
    r, _ = rcen(case, n)
    J, H = core_exact(r)
    rel = t['m1_e']/J - 1.0
    bands = [('r<0.9', r < 0.9), ('0.9-1.1', (r >= 0.9) & (r < 1.1)), ('1.1-3', (r >= 1.1) & (r < 3)),
             ('r>3', r >= 3)]
    fr = t['m1_f1']*r*r
    lim = H[-1]*r[-1]**2   # E units: E <-> J, F/c <-> H (F/cE = H/J)
    st = spl.logstats(d)
    return {'L1': np.mean(np.abs(rel[r > 1.1])), 'txt': 'cyc %d Picard %.2f | E/J_ex-1 ' % (
        t['cycle'], st.get('pic', -1)) + ' '.join('%s %+.2e' % (nm, np.mean(rel[m]))
                                                  for nm, m in bands) +
        ' | F r^2 (r>3) / H_ex r^2 %.4f' % (np.mean(fr[r > 3])/lim)}


def shd_exact(xyz, cen, ca, ns=60):
    """E_ex = Omega_visible / pi of the black-body star r < 1 (flux 1, I = 4 in E units)
    behind an opaque sphere (centre cen, radius ca) from point xyz"""
    r = np.linalg.norm(xyz)
    al = math.asin(min(1.0/r, 1.0))
    mu = 1.0 - (np.arange(ns) + 0.5)/ns*(1.0 - math.cos(al))
    ph = (np.arange(2*ns) + 0.5)/(2*ns)*2.0*math.pi
    M, P = np.meshgrid(mu, ph, indexing='ij')
    st = np.sqrt(1.0 - M*M)
    a = -xyz/r
    t1 = np.cross(a, [0.0, 0.0, 1.0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(a, t1)
    d = (M[..., None]*a + (st*np.cos(P))[..., None]*t1 + (st*np.sin(P))[..., None]*t2)
    d = d.reshape(-1, 3)
    dom = 2.0*math.pi*(1.0 - math.cos(al))/d.shape[0]
    # blocked if the ray hits the clump before the star (clump outside the star)
    oc = xyz - cen
    b = d @ oc
    disc = b*b - (oc @ oc - ca*ca)
    hit = (disc > 0) & (-b - np.sqrt(np.maximum(disc, 0)) > 0)
    return np.sum(~hit)*dom/math.pi, 2.0*(1.0 - math.cos(al))


def ev_shd(case, arm, n, radii=(4.0, 5.0, 5.8)):
    d = run_dir(case, arm, n)
    fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if '.m1.' in f)
    b = spl.binall(os.path.join(d, 'bin', fs[-1]))
    E = b['m1_e']
    n3, n2, n1 = E.shape
    c = spt.case_cfg(case)
    rf = spl.rfaces(c['r'][0], c['r'][1], n1, c['poly'])
    rc = spl.rcent(rf)
    tc = c['th'][0] + (np.arange(n2) + 0.5)*(c['th'][1] - c['th'][0])/n2
    pc = c['ph'][0] + (np.arange(n3) + 0.5)*(c['ph'][1] - c['ph'][0])/n3
    pr = c['prob']
    cr, ct, cp, ca = (float(pr['clump_r']), float(pr['clump_th']), float(pr['clump_ph']),
                      float(pr['clump_a']))
    cen = cr*np.array([math.sin(ct)*math.cos(cp), math.sin(ct)*math.sin(cp), math.cos(ct)])
    out = []
    for rr in radii:
        i = int(np.argmin(np.abs(rc - rr)))
        r = rc[i]
        cls = {'umbra': [], 'penumbra': [], 'lit': []}
        for k in range(n3):
            for j in range(n2):
                u = np.array([math.sin(tc[j])*math.cos(pc[k]), math.sin(tc[j])*math.sin(pc[k]),
                              math.cos(tc[j])])
                if math.acos(min(u @ (cen/cr), 1.0)) > 0.3:
                    continue
                ex, e0 = shd_exact(r*u, cen, ca)
                f = ex/e0
                key = 'umbra' if f < 0.02 else ('lit' if f > 0.98 else 'penumbra')
                cls[key].append((E[k, j, i], ex, e0))
        txt = 'r %.2f' % r
        for key in ('umbra', 'penumbra', 'lit'):
            v = np.array(cls[key])
            if len(v) == 0:
                continue
            if key == 'umbra':
                txt += ' | umbra (%d) <E>/E_free %.3f max %.3f' % (len(v), np.mean(v[:, 0]/v[:, 2]),
                                                                  np.max(v[:, 0]/v[:, 2]))
            else:
                txt += ' | %s (%d) E/E_ex-1 mean %+.3f L1 %.3f' % (
                    key, len(v), np.mean(v[:, 0]/v[:, 1] - 1), np.mean(np.abs(v[:, 0]/v[:, 1] - 1)))
        out.append(txt)
    st = spl.logstats(d)
    return {'L1': 0.0, 'txt': 'Picard %.2f nonconv %g\n    ' % (st.get('pic', -1), st.get('nonconv', -1))
            + '\n    '.join(out)}


def ev_beam(case, arm, n, radii=(1.5, 2.3, 3.0, 4.0, 4.8)):
    """T8: beams injected at r = 1 (vet_gd_beam): exact E = I_b w inside the straight tube
    from the patch along the chosen set direction n (parsed from the log), F = c E n"""
    import re
    d = run_dir(case, arm, n)
    beams = []
    for ln in open(os.path.join(d, 'log')):
        m = re.search(r'vet_gd_beam (\d): direction \d+ n = \(([^)]*)\) w = (\S+) .*I_b = (\S+), '
                      r'patch theta \[(\S+), (\S+)\] phi \[(\S+), (\S+)\].*luminosity \(wedge\) (\S+)', ln)
        if m and int(m.group(1)) not in [b[0] for b in beams]:
            nv = np.array([float(x) for x in m.group(2).split(',')])
            beams.append((int(m.group(1)), nv, float(m.group(3))*float(m.group(4)),
                          [float(m.group(k)) for k in range(5, 9)], float(m.group(9))))
    fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if '.m1.' in f)
    b = spl.binall(os.path.join(d, 'bin', fs[-1]))
    E, F1 = b['m1_e'], b['m1_f1']
    n3, n2, n1 = E.shape
    c = spt.case_cfg(case)
    rf = spl.rfaces(c['r'][0], c['r'][1], n1, c['poly'])
    rc = spl.rcent(rf)
    tf = np.linspace(c['th'][0], c['th'][1], n2 + 1)
    tc = 0.5*(tf[1:] + tf[:-1])
    pc = c['ph'][0] + (np.arange(n3) + 0.5)*(c['ph'][1] - c['ph'][0])/n3
    dom = (np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :]*(c['ph'][1] - c['ph'][0])/n3
    P3, T3 = np.meshgrid(pc, tc, indexing='ij')
    U = np.stack([np.sin(T3)*np.cos(P3), np.sin(T3)*np.sin(P3), np.cos(T3)], -1)
    ltot = sum(bb[4] for bb in beams)
    out = []
    for rr in radii:
        i = int(np.argmin(np.abs(rc - rr)))
        x = rc[i]*U
        ex = np.zeros((n3, n2))
        for _, nv, eb, (t0, t1, p0, p1), _ in beams:
            bq = x @ nv
            disc = bq*bq - (np.sum(x*x, -1) - 1.0)
            sb = bq - np.sqrt(np.maximum(disc, 0.0))     # first r = 1 crossing going back
            p = x - sb[..., None]*nv
            pr = np.linalg.norm(p, axis=-1)
            th = np.arccos(np.clip(p[..., 2]/pr, -1, 1))
            ph = np.arctan2(p[..., 1], p[..., 0])
            ok = (disc > 0) & (sb > 0) & (th >= t0) & (th <= t1) & (ph >= p0) & (ph <= p1)
            ex += np.where(ok, eb, 0.0)
        e = E[:, :, i]
        inb = ex > 0
        lr = np.sum(F1[:, :, i]*dom)*rc[i]**2
        txt = 'r %.2f L/L_inj %.3f' % (rc[i], lr/ltot)
        if inb.any():
            txt += ' | tube (%d cells) E/E_ex median %.3f | E outside tube / total %.3f' % (
                inb.sum(), np.median(e[inb]/ex[inb]), np.sum(e[~inb]*dom[0, None] if False else
                                                              (e*dom)[~inb])/np.sum(e*dom))
        out.append(txt)
    st = spl.logstats(d)
    return {'L1': 0.0, 'txt': 'Picard %.2f nonconv %g, %d beam(s)\n    ' % (
        st.get('pic', -1), st.get('nonconv', -1), len(beams)) + '\n    '.join(out)}


def main():
    mode, case, arm = sys.argv[1:4]
    levs = [int(x) for x in sys.argv[4:]]
    fn = {'thin': ev_thin, 'dif': ev_dif, 'relax': ev_relax, 'core': ev_core,
          'shd': ev_shd, 'beam': ev_beam}[mode]
    res = []
    for n in levs:
        try:
            o = fn(case, arm, n)
        except Exception as ex:  # noqa: BLE001
            print('  %s %s %d: %s' % (case, arm, n, ex))
            continue
        res.append(o)
        print('  %s %s %5d %s' % (case, arm, n, o['txt']))
    if len(res) > 1:
        e = [o['L1'] for o in res]
        print('  orders (L1):', spl.fmto(spl.order(e)))


if __name__ == '__main__':
    main()
