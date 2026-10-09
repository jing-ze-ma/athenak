#!/usr/bin/env python3
"""Steady-cycle kernel cost breakdown from Kokkos Tools simple-kernel-timer (prof batch1).

usage: python3 prof_group.py <labels.tsv> <long dir> [<short dir> <ncyc long> <ncyc short>]
  dir = an arm dir of prof_cuda.sh holding kt_*.txt (kp_reader output, one per rank).
  With a short run the per-kernel times are differenced (long - short): removes setup,
  IC/table reading and the first cycles; times are per rank (mean over the rank files) and
  are printed per steady cycle when the cycle counts are given.
"""
import glob
import os
import re
import sys

# first match wins: (group, label regex, source-file regex)
RULES = [
    ('half-range pass', r'^m1_(vgd_hr|impl_hrq|impl_lfchr)', None),
    ('vet_gd twin', r'^m1_vgd_tw_', None),
    ('vet_gd halo pack/unpack', r'^m1_vgd_(halo|hc_|unpack|seam)', None),
    ('vet_gd sweep (other)', None, r'rad_m1_vetgd'),
    ('vet_col/lat', None, r'rad_m1_(vetcol|vetlat|vet)\.cpp'),
    ('rad implicit halo pack/unpack', r'^m1_(impl_(halo|sct|ncd|hm_)|hi_)', None),
    ('rad mg/line preconditioner', r'^m1_(mg_|gc_|gf_|impl_pcr|impl_thomas|impl_gth)', None),
    ('rad Krylov/BiCGStab vector ops', r'^m1_impl_(bcg|pipe)', None),
    ('rad Krylov other', None, r'rad_m1_(krylov|launch)'),
    ('rad matrix/operator assembly',
     r'^m1_(impl_(op|wb|lag|face|fws|lfc|tbl|od|f1|wk|aphll|ecb|eck|tlim|g0lim|hdir|ktn|src|rhs)'
     r'|vimp_)', None),
    ('rad Picard update/residual/accel',
     r'^m1_(impl_(res|accept|pred|sto|pstore|stop|emax|stb|prein|preout|i0|dgpass|vmmin|thinfrz'
     r'|gms|gcnt|plog)|acc_)', None),
    ('rad other (opacity, explicit flux, time2, coupling, dt)', None, r'rad_m1/'),
    ('halo/MPI pack-unpack (bvals)', None, r'src/bvals/'),
    ('hydro c2p (eos)', None, r'src/eos/'),
    ('hydro fluxes/recon', None, r'src/(hydro/hydro_fluxes|reconstruct/)'),
    ('hydro update/other', None, r'src/hydro/'),
    ('srcterms + pgen (he_star_m1)', None, r'src/(srcterms|pgen)/'),
    ('outputs', None, r'src/outputs/'),
    ('coordinates/mesh/driver/utils', None, r'src/'),
    ('kokkos internal (view init, copies)', r'^Kokkos::', None),
]
RAD = ('half', 'vet', 'rad')
HYD = ('hydro',)

NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def parse(fn):
    """kp_reader text -> ({name: (time, calls)}, total_exec, total_kernels)"""
    k = {}
    tot = ker = None
    name = None
    for line in open(fn, errors='replace'):
        s = line.strip()
        m = re.search(r'Total Execution Time.*?(' + NUM + r')\s*seconds', s)
        if m:
            tot = float(m.group(1))
            continue
        m = re.search(r'Total Time in Kokkos kernels.*?(' + NUM + r')\s*seconds', s)
        if m:
            ker = float(m.group(1))
            continue
        if s.startswith('- '):
            name = s[2:].strip()
            continue
        m = re.match(r'^\((\w+)\)\s+(' + NUM + r')\s+(\d+)', s)
        if m and name is not None:
            t, c = float(m.group(2)), int(m.group(3))
            if m.group(1).lower().startswith('region'):
                name = None
                continue
            a = k.get(name, (0.0, 0))
            k[name] = (a[0] + t, a[1] + c)
            name = None
    return k, tot, ker


def load(d):
    fs = sorted(glob.glob(os.path.join(d, 'kt_*.txt')))
    if not fs:
        sys.exit('no kt_*.txt in %s' % d)
    runs = [parse(f) for f in fs]
    if not any(r[0] for r in runs):
        sys.exit('could not parse %s; first lines:\n%s' % (fs[0], open(fs[0]).read()[:2000]))
    n = len(runs)
    names = set().union(*[r[0] for r in runs])
    mean = {x: (sum(r[0].get(x, (0, 0))[0] for r in runs) / n,
                sum(r[0].get(x, (0, 0))[1] for r in runs) / n) for x in names}
    mx = {x: max(r[0].get(x, (0, 0))[0] for r in runs) for x in names}
    tot = [r[1] for r in runs if r[1] is not None]
    ker = [r[2] for r in runs if r[2] is not None]
    return mean, mx, (sum(tot) / len(tot) if tot else None), \
        (sum(ker) / len(ker) if ker else None), n


def group_of(name, lab):
    f = lab.get(name, '')
    for g, lr, fr in RULES:
        if lr and re.search(lr, name):
            return g
        if fr and f and re.search(fr, f):
            return g
    return 'unmapped (label not a literal in src/)'


def main():
    a = sys.argv[1:]
    if len(a) not in (2, 5):
        sys.exit(__doc__)
    lab = {}
    for line in open(a[0]):
        p = line.rstrip('\n').split('\t')
        if len(p) == 2:
            lab.setdefault(p[0], p[1])
    L, Lmx, Ltot, Lker, nr = load(a[1])
    ncyc = None
    if len(a) == 5:
        S, _, Stot, Sker, _ = load(a[2])
        ncyc = int(a[3]) - int(a[4])
        k = {x: (L[x][0] - S.get(x, (0, 0))[0], L[x][1] - S.get(x, (0, 0))[1]) for x in L}
        tot = Ltot - Stot if (Ltot and Stot) else None
        ker = Lker - Sker if (Lker and Sker) else None
        hdr = 'STEADY = %s minus %s (%d cycles), mean over %d rank files' % (
            a[1], a[2], ncyc, nr)
    else:
        k, tot, ker = L, Ltot, Lker
        hdr = 'WHOLE RUN %s (incl. setup), mean over %d rank files' % (a[1], nr)
    ksum = sum(v[0] for v in k.values())
    ker = ker if ker else ksum
    pc = (lambda t: ' %8.4f s/cyc' % (t / ncyc)) if ncyc else (lambda t: '')
    print('#', hdr)
    print('# wall (Total Execution Time) %s s, in kernels %.3f s (%.1f %%), outside kernels'
          ' (host, MPI waits, launch gaps) %s s' % (
              '%.3f' % tot if tot else 'na', ker, 100 * ker / tot if tot else float('nan'),
              '%.3f' % (tot - ker) if tot else 'na'))
    if ncyc and tot:
        print('# s/cycle (tool run, fenced): %.4f, kernels %.4f' % (tot / ncyc, ker / ncyc))
    print('\n## top 25 kernels (time, % of kernel time, % of wall, calls)')
    top = sorted(k.items(), key=lambda kv: -kv[1][0])[:25]
    for x, (t, c) in top:
        print('%-34s %9.3f s %6.2f %% %6.2f %% %8d  %s%s' % (
            x[:34], t, 100 * t / ker, 100 * t / tot if tot else float('nan'), c,
            group_of(x, lab), pc(t)))
    print('\n## groups (time, % of kernel time, % of wall, kernels)')
    gs = {}
    for x, (t, c) in k.items():
        g = group_of(x, lab)
        gs.setdefault(g, [0.0, 0])
        gs[g][0] += t
        gs[g][1] += 1
    for g, (t, n) in sorted(gs.items(), key=lambda kv: -kv[1][0]):
        print('%-58s %9.3f s %6.2f %% %6.2f %% %4d%s' % (
            g, t, 100 * t / ker, 100 * t / tot if tot else float('nan'), n, pc(t)))
    rad = sum(t for g, (t, n) in gs.items() if g.startswith(RAD))
    hyd = sum(t for g, (t, n) in gs.items() if g.startswith(HYD))
    print('\n## radiation (rad_m1 + vet + half-range) %.3f s = %.1f %% of kernel time;'
          ' hydro (fluxes/update/c2p) %.3f s = %.1f %%; rad/(rad+hydro) = %.1f %%' % (
              rad, 100 * rad / ker, hyd, 100 * hyd / ker, 100 * rad / max(rad + hyd, 1e-30)))
    um = sorted(((x, t) for x, (t, c) in k.items()
                 if group_of(x, lab).startswith('unmapped')), key=lambda kv: -kv[1])[:10]
    if um:
        print('\n## top unmapped labels: ' + ', '.join('%s %.3f s' % u for u in um))
    print('\n## rank spread (long run): top 5 kernels max-rank / mean-rank time')
    for x, (t, c) in top[:5]:
        print('%-34s max %.3f mean %.3f' % (x[:34], Lmx[x], L[x][0]))


if __name__ == '__main__':
    main()
