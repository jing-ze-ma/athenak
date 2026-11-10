"""fe_series.py: incremental fe_hires Kp-marg blueshift series for the dhj page (fehires_fast).

  fe_series.py todo [--stride R] [--max N]   -> prints "planet dump tag" lines still to compute
  fe_series.py run LISTFILE [--gpu]          -> fe3d (both cuts, one trace) + ana per line; keeps the
                                                ana_*.npz in series/<planet>[_mx]/, deletes the spectra
  fe_series.py plot [--img DIR]              -> DIR/fe_vs_rot.png, DIR/fehires_vs_seidel.png,
                                                DIR/fe_takeaway.txt, series/fe_series.json; prints paths
Targets per planet: the newest 5 dumps (the vs-Seidel window) + one dump every R rotations
(R = 10 default), newest first.  Finished runs (W121, W76) stop at rot 300.  Results already
in fehires/out (+ out_mx) are reused (same code, same keys).  Blueshift = -(Kp-marg (a) median);
WASP-76b / WASP-189b combine ESPRESSO and MAROON-X ana by the S/N^2 weight as collect.py.
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
import numpy as np

FF = os.path.dirname(os.path.abspath(__file__))
FE = os.path.join(os.path.dirname(FF), 'fehires')
SER = os.path.join(FF, 'series')
J = '/resnick/groups/carnegie_poc/jingze'
EOS = J + '/w121r1_1101/anchor/eos_1x.txt'
# key, name, Teq, Seidel+2026 blueshift, err, ESPRESSO weight (None = one instrument), run dir, final
PL = [('w76', 'WASP-76b', 2229, 6.95, 0.36, 0.33, J + '/r1planets_1101/wasp76b/runs/p1', True),
      ('w121', 'WASP-121b', 2358, 4.75, 0.30, None, J + '/w121r1_1101/runs/p1', True),
      ('w178', 'WASP-178b', 2398, 4.79, 1.19, None, J + '/r1planets_1101/wasp178b/runs/p1', False),
      ('toi1518', 'TOI-1518b', 2491, 3.82, 0.50, None, J + '/r1planets_1101/toi1518b/runs/p1', False),
      ('w189', 'WASP-189b', 2638, 1.05, 1.59, 0.59, J + '/r1planets_1101/wasp189b/runs/p1', False)]
COL = {'w76': '#b5541c', 'w121': '#1f4e79', 'w178': '#7b3f8c', 'toi1518': '#2e7d4f', 'w189': '#a3822a'}
CUTS = ('pc6', 'pc0')


def dumps(run, final):
    out = {}
    for f in glob.glob(run + '/bin/dhj.hydro_w.*.bin'):
        n = int(re.search(r'\.(\d+)\.bin$', f).group(1))
        if final and n > 600:                  # rot <= 300 (dumps every 0.5 rot)
            continue
        out[n] = f
    return out


def have(key, tag):
    return all(os.path.exists('%s/%s/ana_%s_%s.npz' % (SER, key, tag, c)) for c in CUTS)


def adopt(key, tag, mx):
    """reuse fehires/out (+ out_mx) ana files of the same tag (same code and keys)."""
    src = ['%s/out/ana_%s_%s.npz' % (FE, tag, c) for c in CUTS]
    srcm = ['%s/out_mx/ana_%s_%s.npz' % (FE, tag, c) for c in CUTS]
    if not all(os.path.exists(s) for s in src) or (mx and not all(os.path.exists(s) for s in srcm)):
        return False
    for d, ss in ((key, src), (key + '_mx', srcm if mx else [])):
        if not ss:
            continue
        os.makedirs('%s/%s' % (SER, d), exist_ok=True)
        for s in ss:
            shutil.copy2(s, '%s/%s/%s' % (SER, d, os.path.basename(s)))
    return True


def todo(stride, nmax):
    rows = []
    for key, name, teq, vs, es, wE, run, final in PL:
        ds = dumps(run, final)
        if not ds:
            continue
        ns = sorted(ds)
        for f in glob.glob('%s/out/ana_%s_*_pc6.npz' % (FE, key)):    # every result already made
            t = os.path.basename(f)[4:-8]
            if not have(key, t):
                adopt(key, t, wE is not None)
        want = set(ns[-5:]) | {n for n in ns if n > 0 and n % 2 == 0 and (n//2) % stride == 0}
        for n in sorted(want, reverse=True):
            tag = '%s_%05d' % (key, n)
            if have(key, tag) or adopt(key, tag, wE is not None):
                continue
            rows.append((n not in ns[-5:], -n, key, ds[n], tag))
    rows.sort()                                # newest windows first, then the backfill, newest first
    rows = [r[2:] for r in rows]
    return rows[:nmax] if nmax else rows


def run(listfile, gpu):
    py = sys.executable
    code = os.path.join(FF, 'code')
    tmp = os.path.join(SER, 'tmp')
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for ln in open(listfile):
        if not ln.strip():
            continue
        key, dump, tag = ln.split()
        t0 = time.time()
        cmd = [py, 'fe3d.py', '--planet', 'planets/%s.json' % key, '--eos', EOS, '--out', tmp, dump,
               'tag=' + tag, 'pcuts=1e-6,0'] + (['gpu=1'] if gpu else [])
        with open('%s/fe3d_%s.log' % (tmp, tag), 'w') as lg:
            if subprocess.call(cmd, cwd=code, env=env, stdout=lg, stderr=subprocess.STDOUT):
                print('FAIL fe3d', tag, flush=True)
                continue
        t1 = time.time()
        procs = []
        for c in CUTS:
            for pj, od in ((key, key), (key + '_mx', key + '_mx')):
                if not os.path.exists('%s/planets/%s.json' % (code, pj)):
                    continue
                o = '%s/%s' % (tmp, od)
                os.makedirs(o, exist_ok=True)
                if od != key:
                    s = '%s/spec_%s_%s.npz' % (o, tag, c)
                    if not os.path.exists(s):
                        os.symlink('%s/spec_%s_%s.npz' % (tmp, tag, c), s)
                procs.append((od, c, subprocess.Popen(
                    [py, 'ana.py', '--planet', 'planets/%s.json' % pj, '--out', o if od != key else tmp,
                     '%s_%s' % (tag, c)], cwd=code, env=env, stdout=subprocess.DEVNULL,
                    stderr=open('%s/ana_%s_%s_%s.err' % (tmp, od, tag, c), 'w'))))
        ok = all(p.wait() == 0 for _, _, p in procs)
        for f in glob.glob('%s/ana_*_%s_*.err' % (tmp, tag)):
            if os.path.getsize(f) == 0:
                os.remove(f)
        for od, c, _ in procs:
            src = '%s/ana_%s_%s.npz' % (tmp if od == key else '%s/%s' % (tmp, od), tag, c)
            os.makedirs('%s/%s' % (SER, od), exist_ok=True)
            if os.path.exists(src):
                shutil.move(src, '%s/%s/' % (SER, od))
        for c in CUTS:                          # spectra: 200 MB each, not kept
            for s in glob.glob('%s/spec_%s_%s.npz' % (tmp, tag, c)) + glob.glob('%s/*/spec_%s_%s.npz' % (tmp, tag, c)):
                os.remove(s)
        print('%s %s: fe3d %.0f s, ana %.0f s' % ('done' if ok else 'FAIL ana', tag, t1 - t0, time.time() - t1),
              flush=True)


def blue(key, tag, c, wE):
    """-(Kp-marg (a) 16/50/84) (blueshift > 0), E/MX mixed as collect.py; and rot."""
    z = np.load('%s/%s/ana_%s_%s.npz' % (SER, key, tag, c), allow_pickle=True)
    q = -np.array(z['kpm_full'].item()['q'])[::-1]
    m = '%s/%s_mx/ana_%s_%s.npz' % (SER, key, tag, c)
    if wE is not None and os.path.exists(m):
        q = wE*q + (1 - wE)*(-np.array(np.load(m, allow_pickle=True)['kpm_full'].item()['q'])[::-1])
    return float(z['rot']), q


def series():
    S = {}
    for key, name, teq, vs, es, wE, rdir, final in PL:
        tags = sorted({re.sub(r'_pc[06]\.npz$', '', os.path.basename(f))[4:]
                       for f in glob.glob('%s/%s/ana_%s_*_pc6.npz' % (SER, key, key))})
        tags = [t for t in tags if have(key, t)]
        S[key] = dict(name=name, teq=teq, seidel=vs, err=es, final=final, tags=tags,
                      **{c: [] for c in CUTS}, rot=[])
        for t in tags:
            for c in CUTS:
                r, q = blue(key, t, c, wE)
                S[key][c].append([float(x) for x in q])
            S[key]['rot'].append(r)
    return S


def verdict(s):
    """converged / drifting from the pc6 series: least-squares slope over the last 30 rotations
    (>= 3 points); drifting if |slope| x 30 rot > max(0.2 km/s, 2 x sd of the newest 5 dumps)."""
    rot = np.array(s['rot'])
    if len(rot) < 3:
        return 'too few dumps', np.nan
    b = np.array(s['pc6'])[:, 1]
    o = np.argsort(rot)
    rot, b = rot[o], b[o]
    m = rot >= rot[-1] - 30
    if rot[m][-1] - rot[m][0] < 10:            # needs >= 10 rotations of history
        return 'span < 10 rot', np.nan
    sl = np.polyfit(rot[m], b[m], 1)[0]
    tol = max(0.2, 2*np.std(b[-5:]))
    return ('converged' if abs(sl*30) <= tol else 'drifting'), sl*10


def plot(img):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    os.makedirs(img, exist_ok=True)
    S = series()
    json.dump(S, open(SER + '/fe_series.json', 'w'), indent=1)
    lines = []
    # 1) blueshift vs rotation, one panel per planet (own y-range; Seidel band in each)
    act = [k for k in S if S[k]['rot']]
    fig, axs = plt.subplots(1, len(act), figsize=(3.0*len(act), 3.4), squeeze=False)
    for ax, k in zip(axs[0], act):
        s = S[k]
        o = np.argsort(s['rot'])
        r = np.array(s['rot'])[o]
        ax.axhspan(s['seidel'] - s['err'], s['seidel'] + s['err'], color='#999999', alpha=0.25, lw=0,
                   label='Seidel+2026 1 sigma')
        ax.axhline(s['seidel'], color='#555555', lw=1)
        for c, mf in (('pc6', COL[k]), ('pc0', 'white')):
            q = np.array(s[c])[o]
            ax.plot(r, q[:, 1], '-', color=COL[k], lw=1.2, alpha=0.6 if c == 'pc0' else 1)
            ax.plot(r, q[:, 1], 'o', color=COL[k], mfc=mf, ms=4, mew=1.2,
                    label='Fe cut p<1e-6 bar' if c == 'pc6' else 'no cut')
        v, sl = verdict(s)
        ax.set_title('%s: %s' % (s['name'], v), fontsize=9)
        ax.set_xlabel('rotations', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25, lw=0.5)
        lines.append('%s: %.2f km/s at rot %.1f (1e-6 bar cut; no cut %.2f), %s%s; Seidel %.2f +- %.2f'
                     % (s['name'], np.array(s['pc6'])[o][-1, 1], r[-1], np.array(s['pc0'])[o][-1, 1], v,
                        '' if np.isnan(sl) else ' (%+.2f km/s per 10 rot)' % sl, s['seidel'], s['err']))
    axs[0][0].set_ylabel('Fe I blueshift, Kp-marg (km/s)', fontsize=8)
    axs[0][0].legend(fontsize=7, loc='best')
    fig.tight_layout()
    p1 = img + '/fe_vs_rot.png'
    fig.savefig(p1, dpi=130)
    plt.close(fig)
    # 2) vs Seidel: mean over the newest 5 dumps per planet
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    summ = []
    for k in act:
        s = S[k]
        ax.errorbar(s['teq'], s['seidel'], yerr=s['err'], fmt='o', color='#222222', ms=7, capsize=3, zorder=3)
        o = np.argsort(s['rot'])[-5:]
        r = np.array(s['rot'])[o]
        for c, fill, dx in (('pc6', True, 8), ('pc0', False, 16)):
            Q = np.array(s[c])[o]
            med, lo, hi = Q[:, 1].mean(), Q[:, 0].mean(), Q[:, 2].mean()
            ax.errorbar(s['teq'] + dx, med, yerr=[[med - lo], [hi - med]], fmt='D', color=COL[k],
                        mfc=COL[k] if fill else 'white', mew=1.8, ms=7, capsize=2, zorder=4,
                        label=('%s: rot %.1f-%.1f (%d dumps)%s' % (s['name'], r.min(), r.max(), len(r),
                                                                     '' if s['final'] else ' running'))
                        if c == 'pc6' else None)
            summ.append('%s %s rot %.1f-%.1f n=%d: Kp-marg(a) median %.2f [%.2f, %.2f], dump sd %.2f; Seidel %.2f +- %.2f'
                        % (s['name'], c, r.min(), r.max(), len(r), med, lo, hi, Q[:, 1].std(), s['seidel'], s['err']))
    ax.plot([], [], 'o', color='#222222', label='Seidel+2026 (Fe I net blueshift, 1 sigma)')
    ax.plot([], [], 'D', color='#666666', label='model, Fe cut p<1e-6 bar (open: no cut)')
    ax.set_xlabel('T_eq (K)')
    ax.set_ylabel('Fe I net blueshift (km/s)')
    ax.legend(fontsize=7, loc='upper center', ncol=2)
    ax.set_ylim(-2, max(14, ax.get_ylim()[1]))
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    p2 = img + '/fehires_vs_seidel.png'
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    open(img + '/fe_takeaway.txt', 'w').write('\n'.join(lines) + '\n')
    open(img + '/fehires_vs_seidel_summary.txt', 'w').write('\n'.join(summ) + '\n')
    print(p1)
    print(p2)
    print(img + '/fe_takeaway.txt')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('todo', 'run', 'plot'))
    ap.add_argument('listfile', nargs='?')
    ap.add_argument('--stride', type=int, default=10)
    ap.add_argument('--max', type=int, default=0)
    ap.add_argument('--gpu', action='store_true')
    ap.add_argument('--img', default=os.path.join(SER, 'img'))
    a = ap.parse_args()
    if a.cmd == 'todo':
        for r in todo(a.stride, a.max):
            print(*r)
    elif a.cmd == 'run':
        run(a.listfile, a.gpu)
    else:
        plot(a.img)
