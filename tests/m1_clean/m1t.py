"""m1t.py: small helpers for the clean M1 test suite (Python standard library + numpy).

- read_bin(path): an AthenaK .bin dump as {variable: global array (nx3, nx2, nx1)} plus
  'time' and 'cycle' (single-level meshes; MeshBlocks placed by their logical location).
- read_hst(path): a .hst file as (labels, 2-D array).
- make_deck(src, dst, overrides): copy an input deck, setting block/key = value pairs.
- run(...): run the code in a directory, return (rc, wall seconds, log text).
- compare_dirs(a, b): largest absolute difference over every .bin and .hst file pair.
"""
import glob
import os
import re
import subprocess
import time

import numpy as np


def read_bin(path):
    """Read one AthenaK binary dump (format 1.1 / 1.2) into global arrays."""
    with open(path, 'rb') as fp:
        size = os.fstat(fp.fileno()).st_size
        if fp.readline().split()[0] != b'Athena':
            raise ValueError('%s: not an AthenaK binary dump' % path)
        npar = int(fp.readline().split(b'=')[-1])
        par = {}
        for _ in range(npar - 1):
            k, v = [x.strip() for x in fp.readline().decode().split('=')]
            par[k] = v
        coord_arrays = int(par.get('coordinate arrays', '0')) == 1
        lsz, vsz = int(par['size of location']), int(par['size of variable'])
        fp.readline()   # number of variables
        names = [v.decode() for v in fp.readline().split()[1:]]
        deck = parse_deck_text(fp.read(int(fp.readline().split(b'=')[-1])).decode())
        nx = [int(deck['mesh']['nx%d' % d]) for d in (1, 2, 3)]
        mb = [int(deck['meshblock'].get('nx%d' % d, nx[d - 1])) for d in (1, 2, 3)]
        vdt = np.float64 if vsz == 8 else np.float32
        out = {n: np.zeros((nx[2], nx[1], nx[0])) for n in names}
        while fp.tell() < size:
            ind = np.frombuffer(fp.read(24), dtype=np.int32).astype(np.int64)
            lx = np.frombuffer(fp.read(16), dtype=np.int32)   # logical location, level
            fp.read(6*lsz)                                      # MeshBlock extent
            n1, n2, n3 = ind[1] - ind[0] + 1, ind[3] - ind[2] + 1, ind[5] - ind[4] + 1
            if coord_arrays:
                fp.read((2*(n1 + n2 + n3) + 3)*lsz)
            a = np.fromfile(fp, dtype=vdt, count=n1*n2*n3*len(names))
            a = a.reshape(len(names), n3, n2, n1)
            o1, o2, o3 = lx[0]*mb[0], lx[1]*mb[1], lx[2]*mb[2]
            for v, n in enumerate(names):
                out[n][o3:o3 + n3, o2:o2 + n2, o1:o1 + n1] = a[v]
    out['time'] = float(par['time'])
    out['cycle'] = int(par['cycle'])
    out['_names'] = names
    return out


def parse_deck_text(text):
    """{block: {key: value}} of an input deck's text (comments dropped)."""
    d, blk = {}, None
    for ln in text.splitlines():
        ln = ln.split('#')[0].strip()
        m = re.match(r'<(\w+)>', ln)
        if m:
            blk = m.group(1)
            d.setdefault(blk, {})
        elif blk and '=' in ln:
            k, v = ln.split('=', 1)
            d[blk][k.strip()] = v.strip()
    return d


def read_hst(path):
    """(labels, data) of a history file; labels from its '# [n]=name' header line."""
    labels = []
    with open(path) as f:
        for ln in f:
            if ln.startswith('#') and '[' in ln:
                labels = [x.split('=')[1].strip() for x in ln.split('[')[1:]]
    return labels, np.atleast_2d(np.loadtxt(path))


def make_deck(src, dst, overrides=None, drop=()):
    """Copy deck src to dst with block/key = value overrides (added if absent) and the
    blocks named in drop (e.g. 'output3') removed."""
    want = {}
    for bk, v in (overrides or {}).items():
        b, k = bk.split('/', 1)
        want.setdefault(b, {})[k] = str(v)
    out, blk, seen = [], None, set()
    for ln in open(src).read().splitlines():
        m = re.match(r'\s*<(\w+)>', ln)
        if m:
            if blk in want:
                out += ['%s = %s' % kv for kv in want[blk].items() if (blk, kv[0]) not in seen]
            blk = m.group(1)
            if blk in drop:
                continue
            out.append(ln)
            continue
        if blk in drop:
            continue
        m = re.match(r'\s*(\w+)\s*=', ln)
        if m and blk in want and m.group(1) in want[blk]:
            out.append('%s = %s' % (m.group(1), want[blk][m.group(1)]))
            seen.add((blk, m.group(1)))
            continue
        out.append(ln)
    if blk in want:
        out += ['%s = %s' % kv for kv in want[blk].items() if (blk, kv[0]) not in seen]
    have = set(re.findall(r'^\s*<(\w+)>', '\n'.join(out), re.M))
    for b in want:
        if b not in have:
            out += ['<%s>' % b] + ['%s = %s' % kv for kv in want[b].items()]
    with open(dst, 'w') as f:
        f.write('\n'.join(out) + '\n')


def run(exe, rundir, args, ranks=1, threads=1, mpirun='mpirun', timeout=3600):
    """Run exe with args in rundir; returns (returncode, wall seconds, log text)."""
    os.makedirs(rundir, exist_ok=True)
    cmd = list(exe) + list(args)
    if ranks > 1:
        cmd = mpirun.split() + ['-np', str(ranks)] + cmd
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PROC_BIND='false')
    if ranks == 1:
        # a single-rank run started without mpirun inside a Slurm job: the MPI library
        # would try to join the job's PMI and abort, so it runs as a singleton
        # (Open MPI would also bind the singleton to one core: binding off)
        env = {k: v for k, v in env.items() if not k.startswith(('SLURM', 'PMI')) and
               not (k.startswith('OMPI') and not k.startswith('OMPI_MCA'))}
        env.setdefault('OMPI_MCA_hwloc_base_binding_policy', 'none')
    t0 = time.time()
    with open(os.path.join(rundir, 'run.log'), 'w') as log:
        try:
            rc = subprocess.call(cmd, cwd=rundir, stdout=log, stderr=subprocess.STDOUT,
                                 env=env, timeout=timeout)
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    txt = open(os.path.join(rundir, 'run.log'), errors='replace').read()
    return rc, time.time() - t0, txt


def bins(d, var=None):
    """Sorted .bin dumps under run dir d (optionally of one output variable)."""
    pat = '*.%s.*.bin' % var if var else '*.bin'
    return sorted(glob.glob(os.path.join(d, '**', pat), recursive=True))


def last_bin(d, var):
    fs = bins(d, var)
    if not fs:
        raise FileNotFoundError('no %s dump in %s' % (var, d))
    return read_bin(fs[-1])


def hsts(d):
    return sorted(glob.glob(os.path.join(d, '*.hst')))


def glob_rst(d):
    return sorted(glob.glob(os.path.join(d, '**', '*.rst'), recursive=True))


def _bin_diff(f, g):
    """(max |x - y|, max over variables of max |x - y| / max |x|) of two dumps."""
    x, y = read_bin(f), read_bin(g)
    if x['_names'] != y['_names'] or x['cycle'] != y['cycle']:
        return np.inf, np.inf
    d = [(float(np.abs(x[k] - y[k]).max()), float(np.abs(x[k]).max())) for k in x['_names']]
    return max(a for a, _ in d), max(a/b if b > 0 else a for a, b in d)


def compare_dirs(a, b, last_only=False):
    """(number of file pairs, largest |difference|, largest difference relative to the
    variable's (column's) largest |value|) of the .bin and .hst files of run
    directories a and b.  Default: every file of a against the same file of b.
    last_only (a restarted run b): the last dump of each output variable, and the history
    rows of b against the same number of final rows of a.  A missing partner, a different
    cycle or a different shape counts as an infinite difference."""
    n, mx, mr = 0, 0.0, 0.0
    if last_only:
        vars_ = sorted({os.path.basename(f).split('.')[1] for f in bins(a)})
        pairs = []
        for v in vars_:
            fa, fb = bins(a, v), bins(b, v)
            pairs.append((fa[-1], fb[-1] if fb else None))
    else:
        pairs = [(f, os.path.join(b, os.path.relpath(f, a))) for f in bins(a)]
    for f, g in pairs:
        n += 1
        a_, r_ = (np.inf, np.inf) if not g or not os.path.exists(g) else _bin_diff(f, g)
        mx, mr = max(mx, a_), max(mr, r_)
    for f in hsts(a):
        g = os.path.join(b, os.path.basename(f))
        n += 1
        if not os.path.exists(g):
            return n, np.inf, np.inf
        x, y = read_hst(f)[1], read_hst(g)[1]
        if last_only and y.shape[1] == x.shape[1] and len(y) <= len(x):
            x = x[len(x) - len(y):]
        if x.shape != y.shape:
            return n, np.inf, np.inf
        d = np.abs(x - y).max(axis=0)
        s = np.abs(x).max(axis=0)
        mx = max(mx, float(d.max()))
        mr = max(mr, float(np.max(np.where(s > 0, d/np.where(s > 0, s, 1.0), d))))
    return n, mx, mr


def cell_centres(d):
    """x1, x2, x3 cell centres of the uniform mesh of run directory d's in.athinput."""
    p = parse_deck_text(open(os.path.join(d, 'in.athinput')).read())['mesh']
    out = []
    for k in (1, 2, 3):
        lo, hi, n = float(p['x%dmin' % k]), float(p['x%dmax' % k]), int(p['nx%d' % k])
        out.append(lo + (np.arange(n) + 0.5)*(hi - lo)/n)
    return out


def relax_ode(T0, tot, rho, gm1, rate, t, arad=1.0, nstep=20000):
    """Gas temperature at time t of a uniform gas + radiation box relaxing at constant
    total energy tot = rho T/gm1 + E: de/dt = rate (E - arad T^4), e = rho T/gm1
    (classical RK4)."""
    def f(T):
        return rate*((tot - rho*T/gm1) - arad*T**4)*gm1/rho
    T, h = T0, t/nstep
    for _ in range(nstep):
        k1 = f(T)
        k2 = f(T + 0.5*h*k1)
        k3 = f(T + 0.5*h*k2)
        k4 = f(T + h*k3)
        T += h*(k1 + 2*k2 + 2*k3 + k4)/6.0
    return T


def warnings(txt):
    """The number of CLEAN WARNING lines in a log."""
    return txt.count('CLEAN WARNING')
