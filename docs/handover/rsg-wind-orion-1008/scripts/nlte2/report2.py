"""Tables from thermal_res.pkl / atoms_res.pkl  ->  stdout and tables.txt"""
import re
import sys
import pickle
import numpy as np
HERE = '/orion/ptmp/jinma/rsg_wind_1008/nlte2/'
d = pickle.load(open(HERE + 'thermal_res.pkl', 'rb'))
TG = d['TG']
RES = d['res']
OUT = []


def p(s=''):
    OUT.append(s)


def sel(ic=4, kfac=10., tio='mid', defac=1.):
    return {(r['star'], r['x'], r['P']): r for r in RES
            if r['ic'] == ic and r['kfac'] == kfac and r['tio'] == tio and r['defac'] == defac}


def stab(r, warm=1000.):
    return [q for q in r['roots'] if q['stable']]


def fmt_roots(r):
    s = []
    for q in r['roots']:
        s.append('%.0f%s' % (q['T'], '' if q['stable'] else 'u'))
    tail = '' if (r['H'][0] > 0) else ' (<150: H<0 at 150 K)'
    return '/'.join(s) + tail if s else ('none' + tail)


# previous results from nlte/tables.md (table = lowP, P rows)
prev = {}
for ln in open('/orion/ptmp/jinma/rsg_wind_1008/nlte/tables.md'):
    m = re.match(r'\| (\S+) \| P=(\S+) bar \| (\d) \| ([^|]+) \| ([^|]+) \|', ln)
    if m and (m.group(1), float(m.group(3)), float(m.group(2))) not in prev:
        prev[(m.group(1), float(m.group(3)), float(m.group(2)))] = (m.group(4).strip(),
                                                                    m.group(5).strip())
NOM = sel()
KEYS = sorted(NOM, key=lambda k: (k[0], k[1], k[2]))
tflow = {}
for k in KEYS:
    r = NOM[k]
    x = k[1]
    tflow[k] = (x*r['R']/10e5, x*r['R']/30e5)
TSH = {'golden16': 6.2e7, 'm20lgl5.5': 1.1e8}     # R/v_con from nlte/tables.md

p('## A. Nominal: S_H=0.1, q_forb=1e-11, PI ionisation, v=10 km/s, k_vib(H)=10 k(H2), '
  'TiO mid, dE x1')
p('all roots (u = unstable); stable roots: dominant heating / cooling (share of the total), '
  't_th; LTE and old-model (a) roots from nlte/tables.md')
p()
p('| star | r/R | P bar | roots K | stable root: top heat / top cool | t_th s | t_flow 10/30 s | '
  'R_eff | LTE | old (a) |')
p('|---|---|---|---|---|---|---|---|---|---|')
for k in KEYS:
    r = NOM[k]
    rows = []
    for q in stab(r):
        h = sorted(q['heat'].items(), key=lambda x: -x[1])
        c = sorted(q['cool'].items(), key=lambda x: x[1])
        ht = sum(v for _, v in h)
        ct = -sum(v for _, v in c)
        rows.append('%.0f: %s %.0f%% / %s %.0f%%' % (q['T'], h[0][0], 100*h[0][1]/ht, c[0][0],
                                                     -100*c[0][1]/ct) if h and c else
                    '%.0f' % q['T'])
    tth = ', '.join('%.0e' % q['tth'] for q in stab(r))
    Reff = ', '.join('%.1e' % (q['eps_el']/q['eps_ir']) if q['eps_ir'] != 0 else '-'
                     for q in stab(r))
    pv = prev.get(k, ('?', '?'))
    p('| %s | %.0f | %.0e | %s | %s | %s | %.0e/%.0e | %s | %s | %s |' % (
        k[0], k[1], k[2], fmt_roots(r), '; '.join(rows), tth, *tflow[k], Reff, pv[0], pv[1]))
p()
p('t_sh (R/v_con): golden16 6.2e7 s, m20lgl5.5 1.1e8 s')
p()

# ---------------- brackets
GROUPS = {
    'S_H x q_forb (9)': [dict(ic=i) for i in range(9)],
    'k_vib(H)/k(H2) 1,10,100': [dict(kfac=f) for f in (1., 10., 100.)],
    'TiO lo/mid/hi': [dict(tio=t) for t in ('lo', 'mid', 'hi')],
    'dE_rot x0.5,1,2': [dict(defac=f) for f in (0.5, 1., 2.)],
    'ion LTE (S_H .1 q 1e-11; S_H 1 q 1e-10)': [dict(ic=9), dict(ic=12)],
    'v=30 km/s': [dict(ic=10)],
    'UV(<300nm) x0.01': [dict(ic=11)],
}
p('## B. Brackets around the nominal (one group varied at a time): stable roots, '
  'warm = stable root >= 1000 K; "n/N" = variants with a warm stable root')
p()
hdr = '| star | r/R | P bar | ' + ' | '.join(GROUPS) + ' |'
p(hdr)
p('|' + '---|'*(3 + len(GROUPS)))
for k in KEYS:
    cells = []
    for g, vs in GROUPS.items():
        warm, cold, nw = [], [], 0
        for v in vs:
            a = dict(ic=4, kfac=10., tio='mid', defac=1.)
            a.update(v)
            r = sel(**a)[k]
            st = [q['T'] for q in stab(r)]
            w = [t for t in st if t >= 1000.]
            c = [t for t in st if t < 1000.]
            if not c and r['H'][0] < 0:
                c = [149.]
            warm += w
            cold += c
            nw += bool(w)
        cs = []
        if warm:
            cs.append('w %.0f-%.0f' % (min(warm), max(warm)))
        if cold:
            cs.append('c %s-%.0f' % ('<150' if min(cold) < 150 else '%.0f' % min(cold),
                                     max(cold)) if max(cold) >= 150 else 'c <150')
        cs.append('%d/%d' % (nw, len(vs)))
        cells.append(', '.join(cs))
    p('| %s | %.0f | %.0e | ' % k + ' | '.join(cells) + ' |')
p()
# all-variant census
tot, nw = 0, 0
for r in RES:
    tot += 1
    nw += any(q['stable'] and q['T'] >= 1000 for q in r['roots'])
p('census over all %d variants x conditions: warm stable root present in %d' % (tot, nw))
p()

# ---------------- cooling comparison
sys.path.insert(0, HERE)
import thermal as TH   # noqa: E402
kg = TH.grid_k()
p('## C. Molecular cooling (no radiation, Lambda = n_c n_X L_NK93; per volume, erg cm^-3 s^-1) '
  'vs LTE band emission 4 pi rho sum kl_b B_b(T) (bands 4.4-324 um; all bands), r = 3 R '
  '(golden16, v = 10 km/s, k_vib(H) = 10)')
p()
p('| P bar | T K | n_H2 | n_H | H2O rot | H2O vib | CO rot | CO vib | H2 r+v | SiO r+v | NK sum | '
  'with stellar field (net) | LTE 4.4-324um | LTE all bands |')
p('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
import nlte2 as S1   # noqa: E402
R = NOM[('golden16', 3., 1e-10)]['R']
vr = 10e5/(3*R)
for P in (1e-11, 1e-10, 1e-9):
    iP = list(S1.PBAR).index(P)
    for T0 in (500., 1000., 1500., 2000.):
        it = int(np.argmin(abs(TG - T0)))
        T = TG[it]
        mo0, _, c = TH.molecules(T, P, 1e-30, 4106., vr, 10., 1.)
        mo, _, _ = TH.molecules(T, P, S1.dilution(3.), 4106., vr, 10., 1.)
        Bt = TH.Bb(T)
        rho, kl = kg['rho'][iP][it], kg['kl'][iP][it]
        lf = 4*np.pi*rho*(kl[:3]*Bt[:3]).sum()
        la = 4*np.pi*rho*(kl*Bt).sum()
        p('| %.0e | %.0f | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | '
          '%.1e | %.1e |' % (P, T, c['H2'], c['H'], -mo0['H2O rot'], -mo0['H2O vib'],
                             -mo0['CO rot'], -mo0['CO vib'], -mo0['H2 rot'] - mo0['H2 vib'],
                             -mo0['SiO rot'] - mo0['SiO vib'], -sum(mo0.values()),
                             -sum(mo.values()), lf, la))
p()

# ---------------- eps_Fe
A1 = pickle.load(open(HERE + 'atoms_res.pkl', 'rb'))
p('## D. Fe I / Fe II non-LTE heating per atom h and effective eps (golden16): eps_thin = '
  'h / LTE-thin line heating of the same atom, eps_sob = h / LTE heating x Sobolev beta, '
  'f_abs = h / absorbed (Sobolev) power; x(Fe I) = PI neutral fraction')
p()
p('| q_forb | r/R | P bar | T | n_e | x(Fe I) | Fe I h | Fe I eps_thin | Fe I eps_sob | '
  'Fe II h | Fe II eps_thin | Fe II eps_sob | Fe II f_abs |')
p('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in A1['res']:
    if r['star'] != 'golden16' or r['ic'] not in (3, 4, 5) or r['x'] not in (2., 4.):
        continue
    if r['P'] not in (1e-11, 1e-9):
        continue
    for T0 in (1000., 2000., 3000.):
        q = r['pts'][int(np.argmin(abs(A1['TG'] - T0)))]
        a = q['diag'].get('Fe I', (np.nan,)*5)
        b = q['diag'].get('Fe II', (np.nan,)*5)
        p('| %.0e | %.0f | %.0e | %.0f | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | '
          '%.1e |' % (r['cfg'][1], r['x'], r['P'], q['T'], q['ne'], q['x0']['Fe'], a[0],
                      a[0]/a[1], a[0]/a[2], b[0], b[0]/b[1], b[0]/b[2], b[0]/b[3]))
p()
open(HERE + 'tables.txt', 'w').write('\n'.join(OUT) + '\n')
print('\n'.join(OUT))
