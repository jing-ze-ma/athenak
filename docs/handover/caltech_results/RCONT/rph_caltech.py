"""R_ph and T_eff of hydro_w + vet dumps: restart arms (R1, R25a, R50 rel-f; C50 running binary) against the production dumps.
usage: rph.py STAR   (read only; writes out/STAR.json)
per dump (float32 bins): tau per column from the top face to the cell centre (half cell), kappa_R(rho, T_gas) of the run's Rosseland table;
 RA   = radius where the solid-angle mean tau = 2/3 (linear in log tau between cell centres); RA1 the same for tau = 1 (bsg page definition)
 TA   = (Fcom(RA)/sigma)^(1/4), Fcom = <F1> - 4/3 <E v_r> interpolated in r to RA (shell means)
 RM, TM = medians over the columns of each column's tau = 2/3 radius and of (Fcom_col/sigma)^(1/4) (the wedge page's definition)
"""
import glob, json, os, sys, re
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
J = '/resnick/groups/carnegie_poc/jingze/'
OP2 = J + 'fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
S = dict(m35zams=('prod_c1_1124', 'm35zams', 0.6220, OP2), m35mms=('prod_c1_1124', 'm35mms', 0.6215, OP2),
         bsg=('prod_c1_1124', 'bsg3d', 0.62, J + 'fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt'),
         t9l62=('prod_c4_1124', 't9l62', 0.674, OP2), t9mlt=('prod_mlt_1125', 't9l62', 0.674, OP2))
SIG = 5.670374419e-5
star = sys.argv[1]; pd, base, MU, OP = S[star]; pd = J + pd + '/' + star
NB, NA = int(os.environ.get('NBEF', '20')), int(os.environ.get('NAFT', '2'))
g = [ln for ln in open(OP) if re.match(r'#\s*\d', ln)][0][1:].split()
nT, nD = int(g[0]), int(g[1]); lT0, dlT, lD0, dlD = map(float, g[2:6])
TAB = np.loadtxt(OP, comments='#').reshape(nT, nD)
def kap(rho, T):
    x = np.clip((np.log10(np.maximum(T, 1.0)) - lT0) / dlT, 0, nT - 1 - 1e-9); y = np.clip((np.log10(np.maximum(rho, 1e-30)) - lD0) / dlD, 0, nD - 1 - 1e-9)
    i = x.astype(int); j = y.astype(int); fx = x - i; fy = y - j
    return 10.0 ** ((1 - fx) * (1 - fy) * TAB[i, j] + fx * (1 - fy) * TAB[i + 1, j] + (1 - fx) * fy * TAB[i, j + 1] + fx * fy * TAB[i + 1, j + 1])
def hdr(f):
    h = open(f, 'rb').read(400).decode('latin-1'); return float(re.search(r'time=\s*([-+0-9.eE]+)', h).group(1)), int(re.search(r'cycle=\s*(\d+)', h).group(1))
def asm(d, v):
    lg = np.asarray(d['mb_logical']); b = np.asarray(d['mb_data'][v][0]); n3, n2, n1 = b.shape
    assert int(lg[:, 0].max()) == 0
    P = np.empty((n3 * (lg[:, 2].max() + 1), n2 * (lg[:, 1].max() + 1), n1), np.float64 if np.asarray(b).dtype == np.float64 else np.float32)
    for m in range(len(lg)):
        a, e = lg[m, 1], lg[m, 2]; P[e * n3:(e + 1) * n3, a * n2:(a + 1) * n2] = d['mb_data'][v][m]
    return P
_rows = [l.split() for l in open(pd + '/' + base + '.x1grid.txt') if not l.startswith('#')]; _n = int(re.search(r'is = 3 ie = (\d+)', open(pd + '/' + base + '.x1grid.txt').readline()).group(1)) - 2; gr = np.array([x[:4] for x in _rows if len(x) == 4], float); rc = gr[3:3 + _n, 2]; rf = gr[3:4 + _n, 1]
dr = np.diff(rf)
c = MU * 1.66053906660e-24 / 1.380649e-16
def one(a):
    tag, f = a; cf = f + ''; t, cy = hdr(f)
    h = bc.read_binary(f); ext = (h['x2min'], h['x2max'], h['x3min'], h['x3max'])
    rho = asm(h, 'dens'); T = asm(h, 'eint'); vr = asm(h, 'velx'); dt = str(rho.dtype); del h
    m = bc.read_binary(f.replace('.hydro_w.', '.vet.')); assert abs(m['time'] - t) < 1e-5 * max(1, t)
    E = asm(m, 'vet_e'); F = asm(m, 'vet_f1'); del m
    N3, N2, n1 = rho.shape
    assert n1 == len(rc), (n1, len(rc))
    th = np.linspace(ext[0], ext[1], N2 + 1); W2 = np.cos(th[:-1]) - np.cos(th[1:]); W2 = W2 / W2.sum() / N3
    tm = np.zeros(n1); Fm = np.zeros(n1); Rc = np.empty((N3, N2)); Fc = np.empty((N3, N2))
    for k in range(N3):
        ro = rho[k].astype(float); TT = (2. / 3.) * T[k].astype(float) / ro * c
        d = kap(ro, TT) * ro * dr; tau = np.cumsum(d[:, ::-1], axis=1)[:, ::-1] - 0.5 * d
        fcom = F[k].astype(float) - 4. / 3. * E[k].astype(float) * vr[k].astype(float)
        tm += W2 @ tau; Fm += W2 @ fcom
        th_ = tau >= 2. / 3.; n = np.argmax(th_[:, ::-1], axis=1); ib = np.clip(n1 - 1 - n, 0, n1 - 2)
        j = np.arange(N2); ta, tb = tau[j, ib], tau[j, ib + 1]; fr = np.clip((ta - 2. / 3.) / np.maximum(ta - tb, 1e-30), 0, 1)
        Rc[k] = np.where(th_.any(1), rc[ib] + (rc[ib + 1] - rc[ib]) * fr, np.nan)
        fa, fb = fcom[j, ib], fcom[j, ib + 1]; Fc[k] = fa + (fb - fa) * fr
    def rad(x):
        i = np.where(tm >= x)[0][-1]; a, b = np.log(tm[i]), np.log(tm[i + 1]); return rc[i] + (rc[i + 1] - rc[i]) * (a - np.log(x)) / (a - b)
    RA = rad(2. / 3.); RA1 = rad(1.0); FA = np.interp(RA, rc, Fm)
    ok = np.isfinite(Rc) & (Fc > 0)
    o = dict(tag=tag, f=f, t=t, cyc=cy, dtype=dt, RA=float(RA), RA1=float(RA1), TA=float((FA / SIG) ** 0.25) if FA > 0 else None, LA=float(FA * 4 * np.pi * RA ** 2),
             RM=float(np.nanmedian(Rc)), RMmean=float(np.sum(W2 * Rc)), TM=float(np.median((Fc[ok] / SIG) ** 0.25)), bad=int((~ok).sum()))
    print(star, tag, os.path.basename(f), 'c%d t %.6e' % (cy, t), 'RA %.8e RA1 %.8e TA %s RM %.8e TM %.4f bad %d %s' % (o['RA'], o['RA1'], o['TA'], o['RM'], o['TM'], o['bad'], dt), flush=True)
    return o
if __name__ == '__main__':
    P = sorted(glob.glob(pd + '/bin/' + base + '.hydro_w.*.bin')); A = []
    for a in ('R1', 'R25a', 'R50', 'R25b', 'C50'):
        for f in sorted(glob.glob(J + 'b1_1124/rc/%s/%s/bin/*.hydro_w.*.bin' % (star, a))): A.append((a, f))
    tr = hdr([f for a, f in A if a == 'R1'][0])[0]
    tp = {f: hdr(f)[0] for f in P}; bef = [f for f in P if tp[f] < tr][-NB:]; aft = [f for f in P if tp[f] > tr][:NA]
    aft = [f for f in aft if os.path.exists(f.replace('.hydro_w.', '.vet.')) and __import__('time').time() - os.path.getmtime(f.replace('.hydro_w.', '.vet.')) > 120]
    job = [('P', f) for f in bef] + A + [('Q', f) for f in aft]
    with Pool(int(os.environ.get('NPROC', '4'))) as p: out = p.map(one, job, chunksize=1)
    json.dump(out, open(J + 'rcont_1124/out/%s.json' % star, 'w'), indent=0)
