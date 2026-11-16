#!/usr/bin/env python3
"""compare.py: gate metrics of task_1116_cc_gpu1 runs against the Orion CPU reference.

  python3 compare.py [--runs runs] [--ref ref] [GROUP ...]      gate table (stdout)
  python3 compare.py --make-ref [--runs runs] [--ref ref]       (Orion) write ref/

Per run: rc, FATAL/NaN, CLEAN WARNING and DT_RETRY counts, the CLEAN lines (per-cycle
linear iterations 'sum' and Picard 'pic' against the reference, line by line), the
score (nltek L1 against the exact transport solution, hcyl/cb/shd L1 metrics, rw odd-part
l1_rho/l1_E), conservation (rw: dE, dP; sp: total E and L_z from the final/first bins;
rhse: E_rad + gas, vcs_max from the history), and the field difference of the final
dumps: maxrel = max over variables of max|a - b| / max|b| (blocks assembled into the
global array, so it compares any decomposition).
G2: maxrel of each g2_ run against its g1_ run of the SAME machine (gate <= 1e-7), and
against the CPU 4-rank reference.  The sp 4-block run also against its 1-GPU partner.
Gate widths: score <= max(1e-6, 2 x CPU thread spread), CLEAN 'sum'/'pic' within
max(1, CPU thread spread) per line, decomposition <= max(1e-7, 10 x CPU decomposition
maxrel), WB vcs_max <= max(1e-7, 2 x CPU), dE/dP/dLz <= max(1e-13, 10 x |CPU|) (all in
ref/ref.json: 'spread', 'dec', 'cons')."""
import argparse
import glob
import json
import math
import os
import re
import shutil
import sys

import numpy as np

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/ana')
sys.path.insert(0, H + '/rw')
import bin_convert  # noqa: E402

RW_LINE = ('3.33333 10 x1 32 vetf impl 64 arm=CLN ts=hesdirk2 vimp=true ref=vetqs_d '
           'rm_scheme=clean rm_clean_layout=cell amp=%s rm_vet_gd_rotate_dt=0 '
           'rm_vet_x1_ptol=0')
PARTNER = {'g3_rhse_sp3_mb4': 'g3_rhse_sp3'}   # (cubed sphere: block order differs)


def runs_txt():
    out = []
    for ln in open(H + '/runs.txt'):
        if ln.strip() and not ln.startswith('#'):
            g, tag, n = ln.split()
            out.append((g, tag, int(n)))
    return out


def log_counts(d):
    o = dict(rc=None, bad=0, warn=0, dtr=0, wall=float('nan'), clean=[])
    try:
        txt = open(d + '/run.log', errors='ignore').read()
    except OSError:
        return o
    m = re.findall(r'^rc=(\d+)', txt, re.M)
    o['rc'] = int(m[-1]) if m else None
    o['bad'] = len(re.findall(r'### FATAL|(?<![A-Za-z_])-?(nan|NaN)(?![A-Za-z_])|DIVERGED', txt))
    o['warn'] = txt.count('CLEAN WARNING')
    o['dtr'] = len(re.findall(r'^DT_RETRY', txt, re.M))
    m = re.findall(r'^wall (\S+)', txt, re.M)
    o['wall'] = float(m[-1]) if m else float('nan')
    for c, s, it, sm, pic in re.findall(
            r'^CLEAN cycle (\d+) t2s (\d+) it (-?\d+) \(sum (-?\d+)\) pic (-?\d+)', txt, re.M):
        o['clean'].append([int(c), int(s), int(it), int(sm), int(pic)])
    return o


def last_bins(d):
    """the last dump of every bin output variable in d"""
    out = {}
    for f in glob.glob(d + '/**/*.bin', recursive=True):
        m = re.match(r'.*/([^/]+)\.(\d+)\.bin$', f)
        if m and (m.group(1) not in out or int(m.group(2)) > out[m.group(1)][0]):
            out[m.group(1)] = (int(m.group(2)), f)
    return {k: v[1] for k, v in out.items()}


def first_bins(d):
    out = {}
    for f in glob.glob(d + '/**/*.bin', recursive=True):
        m = re.match(r'.*/([^/]+)\.(\d+)\.bin$', f)
        if m and (m.group(1) not in out or int(m.group(2)) < out[m.group(1)][0]):
            out[m.group(1)] = (int(m.group(2)), f)
    return {k: v[1] for k, v in out.items()}


def assemble(f):
    """{var: global array} of a bin file (uniform mesh; cubed sphere: blocks in gid
    order) and the data dict"""
    b = bin_convert.read_binary(f)
    lg = np.asarray(b['mb_logical'])
    out = {}
    for v in b['var_names']:
        blks = [np.asarray(x) for x in b['mb_data'][v]]
        n3, n2, n1 = blks[0].shape
        if len(set(map(tuple, lg[:, :3]))) == len(blks):
            N = lg[:, :3].max(axis=0) + 1
            g = np.zeros((N[2]*n3, N[1]*n2, N[0]*n1))
            for q, a in enumerate(blks):
                i, j, k = lg[q, 0], lg[q, 1], lg[q, 2]
                g[k*n3:(k+1)*n3, j*n2:(j+1)*n2, i*n1:(i+1)*n1] = a
            out[v] = g
        else:
            out[v] = np.concatenate([a.ravel() for a in blks])
    return out, b


VEC = (('velx', 'vely', 'velz'), ('m1_f1', 'm1_f2', 'm1_f3'), ('mom1', 'mom2', 'mom3'))


def scale(b, v):
    """normalisation of variable v of field dict b: max |b_v|; a vector component by the
    largest magnitude of its vector (velocities at least by max sqrt(eint/dens) ~ c_s),
    so that round-off noise in a near-zero component does not count as a difference"""
    s = np.abs(b[v]).max()
    for grp in VEC:
        if v in grp and all(g in b for g in grp):
            s = max(s, np.sqrt(sum(b[g]**2 for g in grp)).max())
            if grp[0] == 'velx' and 'eint' in b and 'dens' in b:
                s = max(s, np.sqrt(np.abs(b['eint'])/b['dens']).max())
    m = re.match(r'(m1_vh|m1_vf)(\d)$', v)
    if m:
        grp = [k for k in b if re.match(m.group(1) + r'\d$', k)]
        s = max(s, np.sqrt(sum(b[g]**2 for g in grp)).max())
    return s


def maxrel_files(fa, fb):
    a, _ = assemble(fa)
    b, _ = assemble(fb)
    r, worst = 0.0, ''
    for v in b:
        if v not in a or a[v].shape != b[v].shape:
            return float('inf'), v + ' (shape)'
        den = scale(b, v)
        d = np.abs(a[v] - b[v]).max()
        q = d/den if den > 0 else d
        if not np.isfinite(q):
            q = float('inf')
        if q > r:
            r, worst = q, v
    return r, worst


def load_tab(f):
    names = [ln for ln in open(f) if ln.startswith('#')][-1][1:].split()
    a = np.loadtxt(f, comments='#', ndmin=2)
    for x in ('k', 'x3v', 'j', 'x2v'):   # 1-D/2-D tabs list coordinates they do not write
        if len(names) > a.shape[1] and x in names:
            names.remove(x)
    cols = {n: a[:, q] for q, n in enumerate(names)}
    key = [cols[x] for x in ('x1v', 'x2v', 'x3v') if x in cols][::-1]
    o = np.lexsort(key) if key else np.arange(a.shape[0])   # any decomposition
    return {n: v[o] for n, v in cols.items() if n not in ('gid', 'i', 'j', 'k')}


def maxrel_dirs(da, db):
    la, lb = last_bins(da), last_bins(db)
    r, worst = 0.0, ''
    for k in lb:
        if k not in la:
            return float('inf'), k + ' (missing)'
        q, v = maxrel_files(la[k], lb[k])
        if q > r:
            r, worst = q, '%s:%s' % (k, v)
    for ft in sorted(glob.glob(db + '/**/*.tab', recursive=True)):
        stem = re.sub(r'\.\d+\.tab$', '', ft)
        if ft != sorted(glob.glob(stem + '.*.tab'))[-1]:
            continue
        fa = os.path.join(da, os.path.relpath(ft, db))
        if not os.path.exists(fa):
            return float('inf'), os.path.basename(ft) + ' (missing)'
        A, B = load_tab(fa), load_tab(ft)
        for v in B:
            if v not in A or A[v].shape != B[v].shape:
                return float('inf'), '%s:%s (shape)' % (os.path.basename(ft), v)
            den = scale(B, v)
            q = np.abs(A[v] - B[v]).max()/(den if den > 0 else 1.0)
            if q > r:
                r, worst = q, '%s:%s' % (os.path.basename(ft), v)
    return r, worst


def par(d, blk, key):
    txt = open(d + '/in.athinput').read()
    m = re.search(r'<%s>(.*?)(\n<|\Z)' % blk, txt, re.S)
    if m:
        k = re.search(r'^\s*%s\s*=\s*(\S+)' % key, m.group(1), re.M)
        if k:
            return k.group(1)
    return None


def sp_totals(d, which):
    """sp: total energy (gas + radiation) and L_z = sum V R sin(theta) (rho v_phi + F_phi
    /c^2) of the first/last hydro_w + m1 dumps"""
    bs = which(d)
    hk = [k for k in bs if k.endswith('hydro_w')]
    mk = [k for k in bs if k.endswith('.m1')]
    if not hk or not mk:
        return None
    c = float(par(d, 'rad_m1', 'c_light') or 1.0)
    h, _ = assemble(bs[hk[0]])
    m, b = assemble(bs[mk[0]])
    lg = np.asarray(b['mb_logical'])
    N = lg[:, :3].max(axis=0) + 1
    x1a, x1b = float(par(d, 'mesh', 'x1min')), float(par(d, 'mesh', 'x1max'))
    x2a, x2b = float(par(d, 'mesh', 'x2min')), float(par(d, 'mesh', 'x2max'))
    x3a, x3b = float(par(d, 'mesh', 'x3min')), float(par(d, 'mesh', 'x3max'))
    n3, n2, n1 = h['dens'].shape
    rf = np.linspace(x1a, x1b, n1 + 1)
    tf = np.linspace(x2a, x2b, n2 + 1)
    pf = np.linspace(x3a, x3b, n3 + 1)
    del N
    vol = ((pf[1:] - pf[:-1])[:, None, None]*(np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :, None]
           * ((rf[1:]**3 - rf[:-1]**3)/3.0)[None, None, :])
    rc = 0.75*(rf[1:]**4 - rf[:-1]**4)/(rf[1:]**3 - rf[:-1]**3)
    st = np.sin(0.5*(tf[1:] + tf[:-1]))
    R = rc[None, None, :]*st[None, :, None]
    rho = h['dens']
    v2 = h['velx']**2 + h['vely']**2 + h['velz']**2
    eg = h['eint'] + 0.5*rho*v2
    et = float(np.sum(vol*(eg + m['m1_e'])))
    lz = float(np.sum(vol*R*(rho*h['velz'] + m['m1_f3']/c**2)))
    # momentum scale: rho c_s + E/c (the WB decks are static: v, F_phi are noise)
    lzs = float(np.sum(vol*R*(np.sqrt(rho*np.abs(h['eint'])) + m['m1_e']/c))) + 1e-300
    return dict(E=et, Lz=lz, Lzs=lzs)


def hst_totals(d):
    """rhse decks: total energy = gas tot-E (hydro hst) + E_rad (user hst), vcs_max"""
    fh = glob.glob(d + '/*.hydro.hst')
    fu = glob.glob(d + '/*.user.hst')
    if not fh or not fu:
        return None
    def load(f):
        lab = [ln for ln in open(f) if ln.startswith('#')][-1]
        names = [re.sub(r'\[\d+\]=', '', x) for x in lab[1:].split()]
        return names, np.loadtxt(f, comments='#', ndmin=2)
    nh, ah = load(fh[0])
    nu, au = load(fu[0])
    e = ah[:, nh.index('tot-E')] + au[:, nu.index('E_rad')]
    return dict(dE=float(e[-1]/e[0] - 1.0), vcs_max=float(au[-1, nu.index('vcs_max')]),
                vcs_maxall=float(au[1:, nu.index('vcs_max')].max()), E_end=float(e[-1]))


def score(tag, d):
    try:
        if 'nltek' in tag:
            import s2
            from nlte_exact import exact
            m = s2.nlte(d)
            t, J = exact(m['eps'], 6000)
            s = (m['tau'] > 1e-3) & (m['tau'] < 1e5)
            e = np.interp(np.log(m['tau'][s]), np.log(t), J)
            return {'L1': float(np.abs(m['E'][s] - e).sum()/e.sum())}
        for key, dk in (('hcyl', 'hcyl_FIN3_g8x32'), ('cb', 'cb_FIN3_r4'),
                        ('shd', 'shd_FIN3_r12g2')):
            if '_' + key + '_' in tag + '_':
                import score as S
                res = S.metric(dk, d)
                txt = res if isinstance(res, str) else str(res[0])
                nums = re.findall(r'(?<![\w.])[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', txt)
                if not isinstance(res, str):     # (text, the deck's L1 metric)
                    return {'metric': float(res[1])}
                return {'v%d' % q: float(x) for q, x in enumerate(nums)}
    except Exception as e:  # noqa: BLE001
        return {'error': str(e)[:100]}
    return {}


def rw_metrics(runs, tagp, tagm):
    try:
        import rwfarm as F
        import run_radwave as rr
        c, bg, refs, rk, w, vec = F.setup(RW_LINE % '1.00e-03')
        _, period, dt, icfl, nx2, nx3 = rr.make_input(c, bg, w, vec)
        rp, rm = runs + '/' + tagp, runs + '/' + tagm
        out = F.odd_analyse(rp, rm, c, bg, w, vec, nx2, nx3)
        cons = F.conservation(rp, nx2, nx3)
        cons['dmom'] /= 1.0e-3*bg.rho0
        return ({'l1_rho': out['l1_rho'], 'l1_E': out['l1_E'], 'w_re': out['w_re'],
                 'w_im': out['w_im']},
                {'dM': cons['dmass'], 'dE': cons['denergy'], 'dP': cons['dmom']})
    except Exception as e:  # noqa: BLE001
        return {'error': str(e)[:160]}, {}


def metrics(runs, tag):
    d = runs + '/' + tag
    o = log_counts(d)
    o['score'] = score(tag, d)
    if tag.endswith('_rwp'):
        o['score'], o['cons'] = rw_metrics(runs, tag, tag[:-1] + 'm')
    if '_sp' in tag or 'pulsew' in tag:
        a, b = sp_totals(d, first_bins), sp_totals(d, last_bins)
        if a and b:
            o['cons'] = {'dE': b['E']/a['E'] - 1.0, 'dLz': (b['Lz'] - a['Lz'])/b['Lzs']}
    h = hst_totals(d)
    if h:
        o.setdefault('cons', {}).update(h)
    return o


def rel(a, b):
    if a is None or b is None:
        return float('nan')
    return abs(a - b)/max(abs(b), 1e-300)


def clean_diff(a, b):
    """(lines a, lines b, max |d sum|, max |d pic|) over the common CLEAN lines"""
    n = min(len(a), len(b))
    ds = max([abs(a[q][3] - b[q][3]) for q in range(n)] or [0])
    dp = max([abs(a[q][4] - b[q][4]) for q in range(n)] or [0])
    return len(a), len(b), ds, dp


def make_ref(args):
    os.makedirs(args.ref + '/dumps', exist_ok=True)
    ref = {}
    for g, tag, n in runs_txt():
        d = args.runs + '/' + tag
        if not os.path.isdir(d):
            continue
        ref[tag] = metrics(args.runs, tag)
        if g == 'BIT':
            continue
        part = PARTNER.get(tag) or (tag.replace('g2_', 'g1_') if g == 'G2' else None)
        if part and os.path.isdir(args.runs + '/' + part):
            ref[tag]['dec'] = maxrel_dirs(d, args.runs + '/' + part)[0]
        o = args.ref + '/dumps/' + tag
        shutil.rmtree(o, ignore_errors=True)
        os.makedirs(o)
        for k, f in last_bins(d).items():
            shutil.copy(f, o)
        if tag.endswith(('_rwp', '_rwm')) or 'nltek' in tag:   # final tabs only
            for ft in glob.glob(d + '/**/*.tab', recursive=True):
                stem = re.sub(r'\.\d+\.tab$', '', ft)
                if ft == sorted(glob.glob(stem + '.*.tab'))[-1]:
                    q = os.path.join(o, os.path.relpath(os.path.dirname(ft), d))
                    os.makedirs(q, exist_ok=True)
                    shutil.copy(ft, q)
        for f in glob.glob(d + '/*.hst'):
            shutil.copy(f, o)
        shutil.copy(d + '/in.athinput', o)
    # CPU reduction-order spread of the 1-rank runs (OpenMP 1 and 16 threads vs the reference
    # at 4): widens the score / CLEAN-iteration gates where the CPU itself moves (nltek)
    for g, tag, n in runs_txt():
        if g != 'G1' or tag not in ref:
            continue
        sp = {'score': 0.0, 'dsum': 0, 'dpic': 0, 'field': 0.0}
        for alt in args.spread:
            da = alt + '/' + tag
            if not os.path.isdir(da):
                continue
            m = metrics(alt, tag)
            s4 = ref[tag]['score']
            for k, v in m['score'].items():
                if isinstance(v, float) and isinstance(s4.get(k), float):
                    sp['score'] = max(sp['score'], rel(v, s4[k]))
            _, _, ds, dp = clean_diff(m['clean'], ref[tag]['clean'])
            sp['dsum'], sp['dpic'] = max(sp['dsum'], ds), max(sp['dpic'], dp)
            sp['field'] = max(sp['field'], maxrel_dirs(da, args.runs + '/' + tag)[0])
        ref[tag]['spread'] = sp
    json.dump(ref, open(args.ref + '/ref.json', 'w'), indent=1)
    print('wrote', args.ref)


def fmt(x):
    return '%.2e' % x if isinstance(x, float) and np.isfinite(x) else str(x)


def compare(args):
    ref = json.load(open(args.ref + '/ref.json'))
    rows, fails = [], []
    sel = args.groups
    for g, tag, n in runs_txt():
        if g in ('BIT', 'G5') or (sel and g not in sel):
            continue
        d = args.runs + '/' + tag
        if not os.path.isdir(d):
            rows.append('| %s | %s | MISSING |' % (g, tag))
            continue
        m = metrics(args.runs, tag)
        r = ref.get(tag, {})
        la, lb, ds, dp = clean_diff(m['clean'], r.get('clean', []))
        mr, mv = maxrel_dirs(d, args.ref + '/dumps/' + tag)
        sc = {k: rel(v, r.get('score', {}).get(k)) for k, v in m['score'].items()
              if isinstance(v, float)}
        sw = max(sc.values()) if sc else float('nan')
        cons = m.get('cons', {})
        line = ('| %s | %s | %s | %d | %d | %d | %d/%d | %d | %d | %s | %s (%s) | %s | %s |' % (
            g, tag, m['rc'], m['bad'], m['warn'], m['dtr'], la, lb, ds, dp, fmt(sw),
            fmt(mr), mv, ' '.join('%s %s' % (k, fmt(v)) for k, v in cons.items()),
            fmt(m['wall'])))
        rows.append(line)
        bad = []
        if m['rc'] != 0 or m['bad']:
            bad.append('rc/FATAL/NaN')
        if m['warn'] or m['dtr']:
            bad.append('CLEAN WARNING/DT_RETRY')
        sp = r.get('spread', {})
        if la != lb or ds > max(1, sp.get('dsum', 0)) or dp > max(1, sp.get('dpic', 0)):
            bad.append('CLEAN iterations')
        stol = max(1e-6, 2.0*sp.get('score', 0.0))
        if g == 'G1' and sc and not (sw <= stol):
            bad.append('score > %.1e' % stol)
        rc_ = r.get('cons', {})
        if 'vcs_maxall' in cons and not (cons['vcs_maxall'] <= max(1e-7, 2.0*rc_.get(
                'vcs_maxall', 0.0))):
            bad.append('WB vcs_max > max(1e-7, 2 x CPU)')
        for k in ('dE', 'dP', 'dLz'):
            if k in cons and k in rc_ and not (abs(cons[k]) <= max(1e-13, 10*abs(rc_[k]))):
                bad.append('%s not round-off (CPU %.1e)' % (k, rc_[k]))
        part = PARTNER.get(tag) or (tag.replace('g2_', 'g1_') if g == 'G2' else None)
        if part and os.path.isdir(args.runs + '/' + part):
            q, qv = maxrel_dirs(d, args.runs + '/' + part)
            rows.append('|  | %s vs %s (same machine) | maxrel %s (%s) |' % (tag, part, fmt(q), qv))
            lim = max(1e-7, 10.0*r.get('dec', 0.0))
            rows[-1] += ' CPU %s |' % fmt(r.get('dec', float('nan')))
            if not (q <= lim):
                bad.append('decomposition maxrel > %.1e' % lim)
        if bad:
            fails.append('%s: %s' % (tag, ', '.join(bad)))
    print('| grp | run | rc | FATAL/NaN | CLEAN WARN | DT_RETRY | CLEAN lines run/ref | '
          'max d sum | max d pic | score maxrel | field maxrel vs CPU (var) | conservation | '
          'wall s |')
    print('|---' * 15 + '|')
    print('\n'.join(rows))
    print('\nGATE FAILURES:' if fails else '\nall gates passed')
    print('\n'.join(fails))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--runs', default=H + '/runs')
    p.add_argument('--ref', default=H + '/ref')
    p.add_argument('--make-ref', action='store_true')
    p.add_argument('--spread', nargs='*', default=[],
                   help='(make-ref) run roots of the same G1 runs at other thread counts')
    p.add_argument('groups', nargs='*')
    a = p.parse_args()
    a.runs, a.ref = os.path.abspath(a.runs), os.path.abspath(a.ref)
    a.spread = [os.path.abspath(s) for s in a.spread]
    (make_ref if a.make_ref else compare)(a)


if __name__ == '__main__':
    main()
