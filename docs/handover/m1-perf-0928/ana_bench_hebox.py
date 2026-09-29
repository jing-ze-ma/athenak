#!/usr/bin/env python3
"""bench-2026-09-29-hebox: timing and implicit-solve counts from run_bench_hebox.sh logs.

usage: python3 ana_bench_hebox.py <bench.log> [<bench.log> ...] [--c0 C0 --c1 C1]

Per log, the "elapsed=<s> cycle=<n> time=<t> dt=<dt>" lines (ndiag = 10) are cut into
consecutive windows inside [c0, c1] (default 400-1200, skipping the start-up transient;
see README).  Prints one line per run:
  ms/cyc      median over the windows of 1000 * d(elapsed)/d(cycle)
  ms/cyc_w    the same over the whole window (elapsed(c1)-elapsed(c0))/(c1-c0)
  wall/sim-s  median over the windows of d(elapsed)/d(time) [s wall per simulated s]
  wall/sim-s_w the same over the whole window
  dt          mean dt over the window
  pic/solve   Picard passes per implicit transport solve
  kry/solve   BiCGStab inner iterations per transport solve (all its Picard passes)
  kry/lin     BiCGStab inner iterations per linear solve
              (these three come from the <rad_m1> summary lines printed at the end of the
              run, so they cover the WHOLE RUN, cycles 0..nlim)
  NC          NON-CONVERGED implicit transport solves
  fb          fallbacks: positivity + vimp positivity + hesdirk2 stage + line_jacobi
  FATAL count, last cycle, job, host.
"""
import re
import sys
import statistics as st

RE = re.compile(r'^elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)')
NUM = r'([0-9.eE+-]+)'
R_TR = re.compile(r'implicit transport: solves=' + NUM + r' Picard iterations mean=' + NUM
                  + r'.*NON-CONVERGED=' + NUM)
R_BI = re.compile(r'bicgstab: outer passes=' + NUM + r' linear solves=' + NUM
                  + r' .*total=' + NUM)
R_LJ = re.compile(r'bicgstab: breakdowns=' + NUM + r' line_jacobi fallbacks=' + NUM)
R_PO = re.compile(r'closure_lag=\S+ positivity fallbacks=' + NUM)
R_VP = re.compile(r'implicit_vimp positivity fallbacks=' + NUM)
R_H2 = re.compile(r'time_scheme=hesdirk2: stage steps=' + NUM + r' backward-Euler steps='
                  + NUM + r' stage fallbacks=' + NUM)


def ana(path, c0, c1):
    hdr = {}
    rows = []
    nfatal = nnan = 0
    solves = pic = nc = kry = nlin = None
    fb = {}
    for ln in open(path, errors='replace'):
        if ln.startswith('BENCH arm='):
            hdr = dict(kv.split('=', 1) for kv in ln.split()[1:] if '=' in kv)
        m = RE.match(ln)
        if m:
            rows.append((float(m.group(1)), int(m.group(2)), float(m.group(3)),
                         float(m.group(4))))
            continue
        if 'FATAL' in ln:
            nfatal += 1
        if re.search(r'\bnan\b', ln, re.I):
            nnan += 1
        m = R_TR.search(ln)
        if m:
            solves, pic, nc = (float(x) for x in m.groups())
        m = R_BI.search(ln)
        if m:
            nlin, kry = float(m.group(2)), float(m.group(3))
        m = R_LJ.search(ln)
        if m:
            fb['lj'] = float(m.group(2))
        m = R_PO.search(ln)
        if m:
            fb['pos'] = float(m.group(1))
        m = R_VP.search(ln)
        if m:
            fb['vpos'] = float(m.group(1))
        m = R_H2.search(ln)
        if m:
            fb['h2'] = float(m.group(3))
    sel = [r for r in rows if c0 <= r[1] <= c1]
    res = dict(file=path, nr=hdr.get('nranks', '?'), host=hdr.get('host', '?'),
               job=hdr.get('job', '?'), last=rows[-1][1] if rows else -1,
               fatal=nfatal, nan=nnan, nwin=0,
               pic=pic, kry=(kry / solves if kry is not None and solves else None),
               kryl=(kry / nlin if kry is not None and nlin else None),
               nc=nc, fb=(sum(fb.values()) if fb else None))
    if len(sel) < 2:
        return res
    ms, ws = [], []
    for a, b in zip(sel[:-1], sel[1:]):
        dn, de, dtm = b[1] - a[1], b[0] - a[0], b[2] - a[2]
        if dn <= 0 or dtm <= 0:
            continue
        ms.append(1e3 * de / dn)
        ws.append(de / dtm)
    a, b = sel[0], sel[-1]
    res.update(nwin=len(ms), c0=a[1], c1=b[1],
               ms=st.median(ms), ms_w=1e3 * (b[0] - a[0]) / (b[1] - a[1]),
               wps=st.median(ws), wps_w=(b[0] - a[0]) / (b[2] - a[2]),
               dt=(b[2] - a[2]) / (b[1] - a[1]))
    return res


def f(x, fmt):
    return 'NA' if x is None else format(x, fmt)


def main():
    args = sys.argv[1:]
    c0, c1 = 400, 1200
    if '--c0' in args:
        i = args.index('--c0')
        c0 = int(args[i + 1])
        del args[i:i + 2]
    if '--c1' in args:
        i = args.index('--c1')
        c1 = int(args[i + 1])
        del args[i:i + 2]
    print('# ranks ms/cyc ms/cyc_w wall/sim-s wall/sim-s_w dt nwin cycles last '
          'pic/solve kry/solve kry/lin NC fb fatal nan job host file')
    for p in args:
        r = ana(p, c0, c1)
        tail = (f"{f(r['pic'], '.3f')} {f(r['kry'], '.2f')} {f(r['kryl'], '.2f')} "
                f"{f(r['nc'], '.0f')} "
                f"{f(r['fb'], '.0f')} {r['fatal']} {r['nan']} {r['job']} {r['host']} {p}")
        if r['nwin'] == 0:
            print(f"{r['nr']:>3} NO-WINDOW last={r['last']} {tail}")
            continue
        print(f"{r['nr']:>3} {r['ms']:8.2f} {r['ms_w']:8.2f} {r['wps']:9.4f} "
              f"{r['wps_w']:9.4f} {r['dt']:9.4e} {r['nwin']:4d} {r['c0']}-{r['c1']} "
              f"{r['last']} {tail}")


if __name__ == '__main__':
    main()
