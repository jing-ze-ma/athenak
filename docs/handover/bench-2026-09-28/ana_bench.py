#!/usr/bin/env python3
"""bench-2026-09-28: timing from run_bench.sh logs.

usage: python3 ana_bench.py <bench.log> [<bench.log> ...] [--c0 1000 --c1 2000]

Per log: the "elapsed=<s> cycle=<n> time=<t> dt=<dt>" lines (ndiag = 8) are cut into
consecutive windows inside [c0, c1] (default cycles 1000-2000; the start-up transient
and the first Newton calls are skipped).  Prints one line per run:
  ms/cyc   median over the windows of 1000 * d(elapsed)/d(cycle)
  ms/cyc_w the same over the whole window (elapsed(c1)-elapsed(c0))/(c1-c0)
  wall/sim-s median over the windows of d(elapsed)/d(time) [s wall per simulated s]
  wall/sim-s_w the same over the whole window
  dt       mean dt over the window (dt differs slightly between machines)
  nwin, last cycle, FATAL and NOT-CONVERGED counts (the latter needs ck_impl_verbose).
"""
import re
import sys
import statistics as st

RE = re.compile(r'^elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)')


def ana(path, c0, c1):
    hdr = {}
    rows = []
    nfatal = nnc = 0
    for ln in open(path, errors='replace'):
        if ln.startswith('BENCH arm='):
            hdr = dict(kv.split('=', 1) for kv in ln.split()[1:] if '=' in kv)
        m = RE.match(ln)
        if m:
            rows.append((float(m.group(1)), int(m.group(2)), float(m.group(3)),
                         float(m.group(4))))
        if 'FATAL' in ln:
            nfatal += 1
        if 'NOT-CONVERGED' in ln:
            nnc += 1
    sel = [r for r in rows if c0 <= r[1] <= c1]
    res = dict(file=path, arm=hdr.get('arm', '?'), nr=hdr.get('nranks', '?'),
               host=hdr.get('host', '?'), job=hdr.get('job', '?'),
               last=rows[-1][1] if rows else -1, fatal=nfatal, nnc=nnc, nwin=0)
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


def main():
    args = sys.argv[1:]
    c0, c1 = 1000, 2000
    if '--c0' in args:
        i = args.index('--c0')
        c0 = int(args[i + 1])
        del args[i:i + 2]
    if '--c1' in args:
        i = args.index('--c1')
        c1 = int(args[i + 1])
        del args[i:i + 2]
    print('# arm ranks ms/cyc ms/cyc_w wall/sim-s wall/sim-s_w dt nwin cycles last '
          'fatal notconv job host file')
    for p in args:
        r = ana(p, c0, c1)
        if r['nwin'] == 0:
            print(f"{r['arm']:>4} {r['nr']:>3} NO-WINDOW last={r['last']} "
                  f"fatal={r['fatal']} {r['job']} {r['host']} {p}")
            continue
        print(f"{r['arm']:>4} {r['nr']:>3} {r['ms']:8.3f} {r['ms_w']:8.3f} "
              f"{r['wps']:9.5f} {r['wps_w']:9.5f} {r['dt']:8.4f} {r['nwin']:4d} "
              f"{r['c0']}-{r['c1']} {r['last']} {r['fatal']} {r['nnc']} "
              f"{r['job']} {r['host']} {p}")


if __name__ == '__main__':
    main()
