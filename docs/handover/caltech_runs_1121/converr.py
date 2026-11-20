#!/usr/bin/env python3
"""converr.py REF RUN [RUN ...]  (orion CONVERR 2026-11-20; numpy only)

Per-cell solver (iteration) error of one step: every RUN's final state vs the REF's (the tight-tolerance arm),
both from the same restart and the same number of cycles.  Each argument is an armc.sh run root D, with
D/a/bin/*.{hydro_w,m1}.NNNNN.bin (final outputs, written with <output2/3>/bin_double = true),
D/a/run.log and D/a/eos_table.txt (<hydro>/eos_table_dump, the run's own general-EOS table).

Prints, per RUN:
  - bin time/cycle check (the arms must end at the same time: a dt_retry in one arm invalidates the pair)
  - |dE|/E, |dT|/T, |dF|/(cE), |drho|/rho, |deint|/eint: mean, median, p90, p99, p99.9, max (+ where), and
    weighted means (volume; energy = E dV for E and F, eint dV for T; mass for T)
  - the same by region (shell flux ratio f_sh = <F_r>/(c<E>) of REF, a tau proxy ~ 1/(3 f)) and by log-r bin
  - last hst rows: max rel diff per column
  - CLEAN pass lines (rad_m1/clean_nl_print): passes per solve and the last changes of each solve
T(rho, eint) is inverted with a C1 bicubic Hermite fit of the dumped log10(e/rho) nodes (derivatives by finite
differences), Newton + bisection to 1e-14 in log10 T; both arms use the same function, so dT is accurate to
a few per cent of itself even where the fit differs from the code's table at ~1e-5.
"""
import sys, os, glob, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import bin_convert as bc

C_LIGHT = 2.99792458e10   # <units> cgs (length/mass/time_cgs = 1 in the RSG/HG decks)


def final_bin(D, var):
    fs = sorted(glob.glob(os.path.join(D, 'a', 'bin', '*.%s.*.bin' % var)), key=lambda f: int(f.split('.')[-2]))
    if not fs:
        sys.exit('no %s bin in %s/a/bin' % (var, D))
    return fs[-1]


def mbkey(d, m):
    return tuple(np.asarray(d['mb_logical'][m]).tolist()) + tuple(np.round(np.asarray(d['mb_geometry'][m]), 12).tolist())


def varsize(f):
    with open(f, 'rb') as fp:
        h = fp.read(400).decode('ascii', 'ignore')
    m = re.search(r'size of variable=(\d+)', h)
    return int(m.group(1)) if m else -1


# ---------------------------------------------------------------- EOS: T(rho, eint) from the dump
class EOS:
    def __init__(self, f):
        hdr = [l for l in open(f) if l.startswith('#')][:6]
        g = hdr[2].split()[1:]
        self.nx, self.ny = int(g[0]), int(g[1])
        self.x0, self.dx, self.y0, self.dy = map(float, g[2:6])
        a = np.loadtxt(f, usecols=(0,))
        F = a.reshape(self.ny, self.nx)   # F[j, i] = log10(e/rho) at x_i = log10 rho, y_j = log10 T
        self.F = F
        # derivatives in index units (C1 Hermite): central inside, one-sided at the edges
        self.Fx = np.gradient(F, axis=1)
        self.Fy = np.gradient(F, axis=0)
        self.Fxy = np.gradient(self.Fx, axis=0)

    def ev(self, x, y):
        """log10(e/rho) and d/dy (per dex of T) at x = log10 rho, y = log10 T (arrays); linear outside"""
        gx = (x - self.x0) / self.dx
        gy = (y - self.y0) / self.dy
        ix = np.clip(np.floor(gx).astype(np.int64), 0, self.nx - 2)
        iy = np.clip(np.floor(gy).astype(np.int64), 0, self.ny - 2)
        u = gx - ix
        v = gy - iy
        uc = np.clip(u, 0.0, 1.0)
        vc = np.clip(v, 0.0, 1.0)
        u2, u3, v2, v3 = uc * uc, uc ** 3, vc * vc, vc ** 3
        hu = (2 * u3 - 3 * u2 + 1, u3 - 2 * u2 + uc, -2 * u3 + 3 * u2, u3 - u2)
        hv = (2 * v3 - 3 * v2 + 1, v3 - 2 * v2 + vc, -2 * v3 + 3 * v2, v3 - v2)
        dv = (6 * v2 - 6 * vc, 3 * v2 - 4 * vc + 1, -6 * v2 + 6 * vc, 3 * v2 - 2 * vc)
        du = (6 * u2 - 6 * uc, 3 * u2 - 4 * uc + 1, -6 * u2 + 6 * uc, 3 * u2 - 2 * uc)
        f = np.zeros_like(x)
        fy = np.zeros_like(x)
        fx = np.zeros_like(x)
        for a in range(2):
            for b in range(2):
                i0, j0 = ix + a, iy + b
                c = (self.F[j0, i0], self.Fx[j0, i0], self.Fy[j0, i0], self.Fxy[j0, i0])
                # corner (a, b): value slot 2a/2b, slope slot 2a+1 / 2b+1
                ha, hb = hu[2 * a], hv[2 * b]
                ha1, hb1 = hu[2 * a + 1], hv[2 * b + 1]
                da, db = du[2 * a], dv[2 * b]
                da1, db1 = du[2 * a + 1], dv[2 * b + 1]
                f += c[0] * ha * hb + c[1] * ha1 * hb + c[2] * ha * hb1 + c[3] * ha1 * hb1
                fy += c[0] * ha * db + c[1] * ha1 * db + c[2] * ha * db1 + c[3] * ha1 * db1
                fx += c[0] * da * hb + c[1] * da1 * hb + c[2] * da * hb1 + c[3] * da1 * hb1
        fy /= self.dy
        fx /= self.dx
        f = f + (u - uc) * self.dx * fx + (v - vc) * self.dy * fy
        return f, fy

    def logT(self, rho, eint, guess=None):
        x = np.log10(rho)
        g = np.log10(np.maximum(eint, 1e-300) / rho)
        lo = np.full_like(x, self.y0 - 1.0)
        hi = np.full_like(x, self.y0 + (self.ny - 1) * self.dy + 1.0)
        nb = 26
        if guess is not None:   # a nearby state (the REF's T): bracket +-0.01 dex, widened where it does not hold
            f1, _d = self.ev(x, guess - 0.01); f2, _d = self.ev(x, guess + 0.01)
            ok = (f1 <= g) & (f2 > g)
            lo = np.where(ok, guess - 0.01, lo); hi = np.where(ok, guess + 0.01, hi)
            nb = 0 if ok.all() else 26
        for _ in range(nb):   # bisection to ~1e-7 dex
            mid = 0.5 * (lo + hi)
            f, _d = self.ev(x, mid)
            up = f > g
            hi = np.where(up, mid, hi)
            lo = np.where(up, lo, mid)
        y = 0.5 * (lo + hi)
        for _ in range(8):    # Newton, kept in the bracket
            f, fy = self.ev(x, y)
            yn = y - (f - g) / np.where(fy > 0, fy, 1e-30)
            y = np.clip(yn, lo, hi)
        return y


# ---------------------------------------------------------------- geometry (cubed sphere / spherical polar)
def cell_volumes(d, m, nx1, nx2, nx3):
    r = np.asarray(d['mb_x1f'][m])
    dr3 = (r[1:] ** 3 - r[:-1] ** 3) / 3.0
    a = np.asarray(d['mb_x2f'][m]); b = np.asarray(d['mb_x3f'][m])
    if 'use_cubed_sphere' in d.get('_flags', ''):
        s = (np.pi / 4) if max(abs(a).max(), abs(b).max()) > 0.9 else 1.0   # logical [-1,1] or angles
        ta, tb = np.tan(a * s), np.tan(b * s)
        Fc = np.arctan(np.outer(tb, ta) / np.sqrt(1.0 + tb[:, None] ** 2 + ta[None, :] ** 2))   # [k, j]
        dom = Fc[1:, 1:] - Fc[1:, :-1] - Fc[:-1, 1:] + Fc[:-1, :-1]                           # [nx3, nx2]
    else:   # spherical polar: x2 = theta, x3 = phi
        dom = np.outer(np.diff(b), -np.diff(np.cos(a)))
    return dom[:, :, None] * dr3[None, None, :], np.asarray(d['mb_x1v'][m]), dom.sum()


def stats(x, wsets):
    x = np.asarray(x)
    ok = np.isfinite(x)
    xs = x[ok]
    q = np.percentile(xs, [50, 90, 99, 99.9]) if xs.size else [np.nan] * 4
    out = 'mean %.2e med %.2e p90 %.2e p99 %.2e p99.9 %.2e max %.2e' % (xs.mean(), q[0], q[1], q[2], q[3], xs.max())
    for name, w in wsets:
        w = np.asarray(w)[ok]
        out += ' | %s-wt %.2e' % (name, (w * xs).sum() / w.sum())
    if (~ok).sum():
        out += ' | nonfinite %d' % (~ok).sum()
    return out


def hst_cmp(R, D):
    out = []
    for h in sorted(glob.glob(os.path.join(R, 'a', '*.hst'))):
        g = os.path.join(D, 'a', os.path.basename(h))
        if not os.path.exists(g):
            continue
        lab = [l for l in open(h) if l.startswith('#')][-1]
        lab = re.findall(r'\[\d+\]=(\S+)', lab)
        a = [l for l in open(h) if not l.startswith('#')]
        b = [l for l in open(g) if not l.startswith('#')]
        if not a or not b:
            continue
        x = np.array(a[-1].split(), float); y = np.array(b[-1].split(), float); n = min(len(x), len(y))
        rel = np.abs(x[:n] - y[:n]) / np.maximum(np.abs(x[:n]) + np.abs(y[:n]), 1e-300) * 2
        top = np.argsort(-rel)[:8]
        out.append('  hst %s last row (time %s vs %s): ' % (os.path.basename(h), a[-1].split()[0], b[-1].split()[0]) +
                   ', '.join('%s %.2e' % (lab[i] if i < len(lab) else 'c%d' % (i + 1), rel[i]) for i in top))
    return out


def passes(D):
    f = os.path.join(D, 'a', 'run.log')
    if not os.path.exists(f):
        return ['  (no run.log)']
    L = open(f, errors='ignore').read().splitlines()
    pat = re.compile(r'CLEAN pass cycle (\d+) t2s (-?\d+) pos (\d+) pass (\d+) it (-?\d+) de (\S+) dv (\S+)')
    solves = []
    for l in L:
        m = pat.search(l)
        if not m:
            continue
        c, s, p, n, it = map(int, m.groups()[:5]); de, dv = float(m.group(6)), float(m.group(7))
        if n == 0 or not solves:
            solves.append([])
        solves[-1].append((c, s, p, n, it, de, dv))
    nret = sum(int(m.group(1)) for l in L for m in [re.search(r'dt_retry: \d+ checkpoints .*?, (\d+) retries', l)] if m)
    cyc = [l.split()[1] for l in L if l.startswith('elapsed=')]
    out = ['  run.log: CW %d FATAL %d dt_retry retries %d; cycle lines %s' % (sum('CLEAN WARNING' in l for l in L),
                                                                  sum('FATAL' in l for l in L), nret, ' '.join(cyc))]
    if solves:
        np_ = [len(s) for s in solves]
        out.append('  CLEAN pass: %d solves, passes/solve %s (total %d)' % (len(solves), np_, sum(np_)))
        for s in solves:
            c, t2, p = s[0][:3]
            seq = ' '.join('%.1e/%.1e(%d)' % (x[5], x[6], x[4]) for x in s[-4:])
            out.append('    cycle %d t2s %d pos %d: %d passes, last de/dv(it): %s' % (c, t2, p, len(s), seq))
    else:
        out.append('  (no CLEAN pass lines: run without rad_m1/clean_nl_print)')
    return out


def load(D):
    h = bc.read_binary(final_bin(D, 'hydro_w'))
    r = bc.read_binary(final_bin(D, 'm1'))
    return h, r


def main():
    R, runs = sys.argv[1], sys.argv[2:]
    eosf = os.path.join(R, 'a', 'eos_table.txt')
    eos = EOS(eosf)
    print('REF %s (EOS %s: %d x %d nodes)' % (R, eosf, eos.nx, eos.ny))
    hR, mR = load(R)
    for f in (final_bin(R, 'hydro_w'), final_bin(R, 'm1')):
        print('  REF bin %s time %.12g size of variable %d' % (os.path.basename(f), hR['time'], varsize(f)))
    hdr = hR['header'] if isinstance(hR['header'], str) else ' '.join(map(str, hR['header']))
    flags = 'use_cubed_sphere' if re.search(r'use_cubed_sphere\s*=\s*(true|1)\b', hdr) else ''
    hR['_flags'] = flags
    print('  mesh: %s, %d MeshBlocks, coordinate arrays %s' % ('cubed sphere' if flags else 'spherical polar', hR['n_mbs'],
                                                             'yes' if 'mb_x1f' in hR else 'NO'))
    hv = hR['var_names']; mv = mR['var_names']
    iD, iE = hv[0], hv[4]
    ordR = sorted(range(hR['n_mbs']), key=lambda m: mbkey(hR, m))
    ordRm = sorted(range(mR['n_mbs']), key=lambda m: mbkey(mR, m))
    # REF per-cell arrays (flattened in a fixed MeshBlock order)
    V, RV, E, F1, Fv, RHO, EI, LT, IR = [], [], [], [], [], [], [], [], []
    omega = 0.0
    for m, mm in zip(ordR, ordRm):
        dv, rv, om = cell_volumes(hR, m, 0, 0, 0)
        omega += om
        V.append(dv.ravel()); RV.append(np.broadcast_to(rv[None, None, :], dv.shape).ravel())
        IR.append(np.broadcast_to(np.arange(dv.shape[2])[None, None, :], dv.shape).ravel())
        e = np.asarray(mR['mb_data']['m1_e'][mm], float)
        E.append(e.ravel()); F1.append(np.asarray(mR['mb_data']['m1_f1'][mm], float).ravel())
        Fv.append(np.stack([np.asarray(mR['mb_data']['m1_f%d' % q][mm], float).ravel() for q in (1, 2, 3)]))
        rho = np.asarray(hR['mb_data'][iD][m], float).ravel(); ei = np.asarray(hR['mb_data'][iE][m], float).ravel()
        RHO.append(rho); EI.append(ei); LT.append(eos.logT(rho, ei))
    V = np.concatenate(V); RV = np.concatenate(RV); IR = np.concatenate(IR); E = np.concatenate(E)
    F1 = np.concatenate(F1); Fv = np.concatenate(Fv, axis=1); RHO = np.concatenate(RHO); EI = np.concatenate(EI)
    LT = np.concatenate(LT)
    nr = IR.max() + 1
    print('  total solid angle / 4pi = %.6f (cell-volume check; 1 = full sphere)' % (omega / (4 * np.pi)))
    # shell flux ratio (REF): f_sh = sum F_r dV / (c sum E dV), per radial index
    fsh = np.bincount(IR, F1 * V, nr) / (C_LIGHT * np.bincount(IR, E * V, nr))
    rsh = np.bincount(IR, RV * V, nr) / np.bincount(IR, V, nr)
    Tsh = 10 ** (np.bincount(IR, LT * V, nr) / np.bincount(IR, V, nr))
    print('  shells: r %.4e .. %.4e cm; f_sh(0) %.2e; r(f_sh>0.03) %.4e, r(f_sh>0.1) %.4e, r(f_sh>0.3) %.4e' % (
        rsh[0], rsh[-1], fsh[0],
        *[rsh[np.argmax(fsh > t)] if (fsh > t).any() else np.nan for t in (0.03, 0.1, 0.3)]))
    REG = np.digitize(fsh, [1e-3, 3e-2, 0.3])[IR]
    regn = ['interior f<1e-3 (tau>~300)', 'envelope 1e-3..0.03', 'photosphere 0.03..0.3 (tau~1-10)', 'wind/thin f>0.3']
    lr = np.log10(rsh); edges = np.linspace(lr[0], lr[-1] + 1e-12, 9)
    RB = np.digitize(np.log10(RV), edges) - 1
    for D in runs:
        print('\n=== RUN %s vs REF' % D)
        hD, mD = load(D)
        tR, tD = hR['time'], hD['time']
        print('  time REF %.15g RUN %.15g (%s); cycle %s vs %s' % (tR, tD, 'SAME' if tR == tD else 'DIFFERENT: pair invalid',
                                                               hR['cycle'], hD['cycle']))
        ordD = sorted(range(hD['n_mbs']), key=lambda m: mbkey(hD, m))
        ordDm = sorted(range(mD['n_mbs']), key=lambda m: mbkey(mD, m))
        if [mbkey(hD, m) for m in ordD] != [mbkey(hR, m) for m in ordR]:
            print('  MeshBlock sets differ: skipped'); continue
        dE, dT, dF, dRho, dEi = [], [], [], [], []
        k = 0
        for m, mm in zip(ordD, ordDm):
            e = np.asarray(mD['mb_data']['m1_e'][mm], float).ravel(); n = e.size; s = slice(k, k + n)
            dE.append(np.abs(e - E[s]) / np.abs(E[s]))
            fv = np.stack([np.asarray(mD['mb_data']['m1_f%d' % q][mm], float).ravel() for q in (1, 2, 3)])
            dF.append(np.sqrt(((fv - Fv[:, s]) ** 2).sum(0)) / (C_LIGHT * np.abs(E[s])))
            rho = np.asarray(hD['mb_data'][iD][m], float).ravel(); ei = np.asarray(hD['mb_data'][iE][m], float).ravel()
            dT.append(np.abs(10 ** (eos.logT(rho, ei, LT[s]) - LT[s]) - 1.0))
            dRho.append(np.abs(rho - RHO[s]) / RHO[s]); dEi.append(np.abs(ei - EI[s]) / EI[s])
            k += n
        dE, dT, dF, dRho, dEi = map(np.concatenate, (dE, dT, dF, dRho, dEi))
        WE, WT, WM = E * V, EI * V, RHO * V
        print('  cells %d; identical E %d, identical T %d' % (dE.size, (dE == 0).sum(), (dT == 0).sum()))
        for nm, x, ws in (('|dE|/E      ', dE, [('vol', V), ('E', WE)]),
                          ('|dT|/T      ', dT, [('vol', V), ('eint', WT), ('mass', WM)]),
                          ('|dF|/(cE)   ', dF, [('vol', V), ('E', WE)]),
                          ('|drho|/rho  ', dRho, [('vol', V), ('mass', WM)]),
                          ('|deint|/eint', dEi, [('vol', V), ('eint', WT)])):
            print('  ALL %s %s' % (nm, stats(x, ws)))
        for nm, x in (('dE/E', dE), ('dT/T', dT), ('dF/cE', dF)):
            i = np.nanargmax(x)
            print('    max %s %.2e at r %.4e (shell %d, f_sh %.2e, T %.3e K, rho %.2e, E %.2e, |F|/cE %.2e)' % (
                nm, x[i], RV[i], IR[i], fsh[IR[i]], 10 ** LT[i], RHO[i], E[i], np.sqrt((Fv[:, i] ** 2).sum()) / (C_LIGHT * E[i])))
        print('  by region (REF shell flux ratio f_sh):')
        for g in range(4):
            sel = REG == g
            if not sel.any():
                continue
            rr = RV[sel]
            print('   [%s] cells %d, r %.3e..%.3e, T %.2e..%.2e K' % (regn[g], sel.sum(), rr.min(), rr.max(),
                                                                 10 ** LT[sel].min(), 10 ** LT[sel].max()))
            print('      dE/E  %s' % stats(dE[sel], [('vol', V[sel]), ('E', WE[sel])]))
            print('      dT/T  %s' % stats(dT[sel], [('vol', V[sel]), ('eint', WT[sel])]))
            print('      dF/cE %s' % stats(dF[sel], [('vol', V[sel]), ('E', WE[sel])]))
        print('  by log r bin (8 bins):')
        for b in range(8):
            sel = RB == b
            if not sel.any():
                continue
            ish = np.unique(IR[sel])
            print('   r %.3e..%.3e (shells %d-%d, f_sh %.1e..%.1e, T_sh %.2e..%.2e): dE/E mean %.2e p99 %.2e max %.2e | '
                  'dT/T mean %.2e p99 %.2e max %.2e | dF/cE mean %.2e max %.2e' % (
                      10 ** edges[b], 10 ** edges[b + 1], ish.min(), ish.max(), fsh[ish].min(), fsh[ish].max(),
                      Tsh[ish].min(), Tsh[ish].max(), dE[sel].mean(), np.percentile(dE[sel], 99), dE[sel].max(),
                      dT[sel].mean(), np.percentile(dT[sel], 99), dT[sel].max(), dF[sel].mean(), dF[sel].max()))
        for l in hst_cmp(R, D):
            print(l)
        for l in passes(D):
            print(l)
        del hD, mD
    print('\n=== REF passes')
    for l in passes(R):
        print(l)


if __name__ == '__main__':
    main()
