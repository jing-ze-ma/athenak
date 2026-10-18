"""Photosphere seam removal (tools/ic_star, 10-16).

(1) Scaffold TAPER: the frozen F_MLT/F of the column is multiplied by a smoothstep in ln tau,
    1 at tau >= tau_hi (30), 0 at tau <= tau_lo (3); F_r = F (1 - F_MLT/F) keeps F_r + F_MLT =
    L exactly.  Why tau 30 -> 3: he_ic_rad = sph1d replaces T by the radiative-equilibrium T
    with weight w = 1 for tau < he_ic_rad_taua = 30 (he_star_m1.cpp ~821), i.e. above tau 30 the
    code's IC is RE with the column's F_r; a scaffold that is still ~0.5 F there (He giant) or
    ends abruptly at the column's photosphere (the hg145/BSG gate failures) is inconsistent
    with it.  Tapering inside the full-RE zone hands the whole flux to radiation by tau 3, the
    Hopf join of the column, smoothly.
(2) Code-consistent atmosphere: the column is passed through the code's OWN startup solve
    (he_ic_rad = sph1d + re_atm + he_ic_rad_hse, the run's Rosseland/Planck tables, tau_a/tau_b
    blend, HsSph1dPass) by running the production binary (CPU build of the same commit) for 0
    cycles with problem/he_ic_rad_dump; the dumped column (+ the re-applied taper) becomes the
    new column; iterate to the fixed point R_ph(column) = R_ph(sph1d) (the code prints
    'R_ph (tau_R 2/3, run kappa) rph0 -> rph1').  Then the run's startup reproduces the column's
    photosphere to within the tolerance (default 0.1 of the radial cell at R_ph).
"""
import os
import re
import subprocess

import numpy as np



def smoothstep_tau(tau, lo, hi):
    x = (np.log(np.maximum(tau, 1e-300)) - np.log(lo))/(np.log(hi) - np.log(lo))
    x = np.clip(x, 0.0, 1.0)
    return x*x*(3.0 - 2.0*x)


def tau_from_top(r, rho, kap, rtop):
    """tau from the mesh top (rows above rtop contribute 0), trapezoid: as the pgen"""
    a = np.where(r < rtop, rho*kap, 0.0)
    tau = np.zeros(r.size)
    tau[:-1] = np.cumsum((0.5*(a[1:] + a[:-1])*np.diff(r))[::-1])[::-1]
    return tau


def taper(r, rho, T, fm0, L, kap_fn, rtop, lo, hi, extend=False):
    """-> (fmlt, F_r, tau): fm0 tapered by the smoothstep in ln tau.  extend: the scaffold is
    first continued outward at its value at tau_hi (for a column whose F_MLT ends at the Hopf
    join tau 3 and must not be removed below it), then tapered to 0 at tau_lo"""
    tau = tau_from_top(r, rho, kap_fn(T, rho), rtop)
    s = smoothstep_tau(tau, lo, hi)
    if extend:
        k = np.where(tau >= hi)[0]
        if k.size:
            fm0 = np.where(tau < hi, fm0[k[-1]], fm0)
    fm = fm0*s
    F = L/(4*np.pi*r**2)
    return fm, F*(1.0 - fm), tau


def write_col(fn, hdr, r, rho, eint, Fr, E, T, fm):
    with open(fn, 'w') as fh:
        for h in hdr:
            fh.write(h if h.endswith('\n') else h + '\n')
        fh.write('# r[cm] rho eint F_r(radiative part) E  | T[K] F_MLT/F (cols 6,7)\n')
        for row in zip(r, rho, eint, Fr, E, T, fm):
            fh.write('%.10e %.10e %.10e %.10e %.10e %.10e %.6e\n' % row)


def read_col(fn):
    d = np.loadtxt(fn, comments='#')
    return d.T[:7]


RPH_RE = re.compile(r'R_ph \(tau_R 2/3, run kappa\) ([0-9.eE+-]+) -> ([0-9.eE+-]+)')


def deck_with(base, keys, out):
    """copy of the deck file with block/key = value set (replaced, or inserted after the
    block header when the key is not there)"""
    lines = open(base).read().splitlines()
    todo = {}
    for kv in keys:
        k, v = kv.split('=', 1)
        b, kk = k.split('/')
        todo.setdefault(b, {})[kk] = v
    res, blk, seen = [], None, set()
    for ln in lines:
        st = ln.strip()
        if st.startswith('<') and st.endswith('>'):
            blk = st[1:-1]
            res.append(ln)
            for kk, v in todo.get(blk, {}).items():
                if not any(x.split('=')[0].strip() == kk for x in _block(lines, blk)):
                    res.append('%s = %s' % (kk, v))
                    seen.add((blk, kk))
            continue
        name = st.split('=')[0].strip() if '=' in st and not st.startswith('#') else None
        if blk in todo and name in todo[blk]:
            res.append('%s = %s' % (name, todo[blk][name]))
            seen.add((blk, name))
        else:
            res.append(ln)
    open(out, 'w').write('\n'.join(res) + '\n')
    return out


def _block(lines, blk):
    out, on = [], False
    for ln in lines:
        st = ln.strip()
        if st.startswith('<') and st.endswith('>'):
            on = (st[1:-1] == blk)
            continue
        if on and '=' in st and not st.startswith('#'):
            out.append(st)
    return out


def code_pass(cfg, colfile, dumpfile, wdir, log):
    """0-cycle run of the binary on the deck with the column; -> rph0, rph1, stdout path"""
    os.makedirs(wdir, exist_ok=True)
    keys = ['problem/he_ic_file=%s' % colfile, 'problem/he_ic_rad_dump=%s' % dumpfile,
            'time/nlim=0', 'problem/he_ic_equil=false', 'rad_m1/ic_rad=none',
            'problem/he_ic_rad_atm_force=continuum'] + list(cfg.get('overrides', []))
    deck = deck_with(cfg['deck'], keys, os.path.join(wdir, 'pass.athinput'))
    cmd = list(cfg.get('launcher', [])) + [cfg['binary'], '-i', deck, '-d', wdir]
    out = os.path.join(wdir, 'run.log')
    env = dict(os.environ, OMP_NUM_THREADS=str(cfg.get('omp', 8)), OMP_PROC_BIND='false')
    with open(out, 'w') as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=wdir,
                            timeout=cfg.get('timeout', 7200)).returncode
    txt = open(out).read()
    m = RPH_RE.search(txt)
    if rc != 0 or m is None or not os.path.exists(dumpfile):
        raise RuntimeError('seam: the code pass failed (rc %d), see %s:\n%s'
                           % (rc, out, txt[-2000:]))
    h = re.search(r'HSE last max \|dln rho\| ([0-9.eE+-]+).*?in (\d+) fine cells', txt)
    log('    code pass: %s; HSE last max |dln rho| %s, Gamma-clamped fine cells %s%s'
        % (out, h.group(1) if h else '?', h.group(2) if h else '?',
           ' -- HSE NOT CONVERGED' if 'did NOT converge' in txt else ''))
    return float(m.group(1)), float(m.group(2)), out


def grid_dr_at(wdir, r0):
    """the radial cell width at r0 from <basename>.x1grid.txt (he_grid_dump)"""
    fs = [f for f in os.listdir(wdir) if f.endswith('.x1grid.txt')]
    if not fs:
        return None
    d = np.array([[float(v) for v in ln.split()[1:4]] for ln in open(os.path.join(wdir, fs[0]))
                  if not ln.startswith('#') and len(ln.split()) >= 4])
    i = int(np.argmin(abs(d[:, 1] - r0)))
    return float(d[i, 2])


def fixed_point(cfg, col_r, col_fm0, L, kap_fn, rtop, hdr, base, log, lo, hi, extend=False):
    """iterate column -> code sph1d dump -> taper -> column to R_ph(col) = R_ph(sph1d).
    col_r etc.: the starting column file 'base' (already tapered); col_fm0: the UNTAPERED
    F_MLT/F as a function of r (callable)."""
    wd = cfg['workdir']
    tol = cfg.get('tol_cells', 0.1)
    hist = []
    cur = base
    for it in range(cfg.get('maxit', 8)):
        dump = os.path.join(wd, 'dump_%d.txt' % it)
        rph0, rph1, out = code_pass(cfg, cur, dump, os.path.join(wd, 'pass_%d' % it), log)
        dr = grid_dr_at(os.path.join(wd, 'pass_%d' % it), rph1) or np.nan
        hist.append((it, rph0, rph1, dr))
        log('  seam pass %d: R_ph(column) %.8e -> R_ph(sph1d) %.8e: dR %.4e cm = %.4f cells '
            '(dr %.4e)' % (it, rph0, rph1, rph1 - rph0, (rph1 - rph0)/dr, dr))
        if abs(rph1 - rph0) < tol*dr and it > 0:
            return cur, hist
        r, rho, eint, Fr, E, T, fm = read_col(dump)
        # the dump spans exactly the fine grid; keep the previous column's rows outside it
        # (the pgen needs strict cover of the mesh + ghost margin)
        b = np.array(read_col(cur))
        lo_, hi_ = b[:, b[0] < r[0]*(1 - 1e-9)], b[:, b[0] > r[-1]*(1 + 1e-9)]
        r, rho, eint, Fr, E, T, fm = [np.concatenate([lo_[j], v, hi_[j]]) for j, v in
                                      enumerate((r, rho, eint, Fr, E, T, fm))]
        fmt, Frt, tau = taper(r, rho, T, col_fm0(r), L, kap_fn, rtop, lo, hi, extend)
        nxt = os.path.join(wd, 'col_%d.txt' % (it + 1))
        write_col(nxt, hdr + ['# seam pass %d: rows = the code dump %s (sph1d re_atm hse), '
                              'taper tau %g..%g re-applied' % (it, dump, lo, hi)],
                  r, rho, eint, Frt, E, T, fmt)
        cur = nxt
    log('  seam: NOT converged in %d passes' % len(hist))
    return cur, hist
