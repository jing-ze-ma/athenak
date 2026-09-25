"""tmr.py <run dir> [c0=10]: implicit_timers breakdown in ms/cycle (cycles c0..end) and
the elapsed wall over the same window (cycle c0 to the last cycle line)."""
import re
import sys
d = sys.argv[1]
c0 = int(sys.argv[2]) if len(sys.argv) > 2 else 10
t = {}
tm = {}
cnt = {}
for ln in open(d + '/run.log'):
    m = re.search(r'elapsed=([0-9.e+-]+) cycle=(\d+)', ln)
    if m:
        t[int(m.group(2))] = float(m.group(1))
    if 'timers(s):' in ln:
        for k, v in re.findall(r'(\w+)=([0-9.e+-]+)', ln.split('timers(s):')[1]):
            tm[k] = float(v)
    if 'timers from cycle' in ln:
        for k, v in re.findall(r'(\w+)=([0-9.e+-]+)', ln):
            cnt[k] = float(v)
cl = max(t)
last = max(int(x) for x in re.findall(r'cycle=(\d+)', open(d + '/run.log').read()))
ncyc = last - c0
c1 = int(sys.argv[3]) if len(sys.argv) > 3 else 50
wall = (t[c1] - t[c0])/(c1 - c0)*1e3
rad = sum(tm.values())*1e3/ncyc
print(f'{d}: {ncyc:.0f} cycles timed; wall {wall:.2f} ms/cycle (cycles {c0}-{c1}, '
      f'includes fences); radiation {rad:.2f} ms/cycle; rest (hydro etc.) '
      f'{wall - rad:.2f}; rad/rest {rad/(wall - rad):.2f}')
for k, v in tm.items():
    print(f'   {k:11s} {v*1e3/ncyc:7.2f} ms/cycle')
s = cnt.get('solves', 0)
print(f'   per cycle: solves {s/ncyc:.2f} passes {cnt.get("passes",0)/ncyc:.2f} '
      f'krylov it {cnt.get("kry_it",0)/ncyc:.1f}; per stage: passes s1 '
      f'{cnt.get("passes_s1",0)/ncyc:.2f} s2 {cnt.get("passes_s2",0)/ncyc:.2f}, '
      f'it s1 {cnt.get("it_s1",0)/ncyc:.1f} s2 {cnt.get("it_s2",0)/ncyc:.1f}; '
      f'krylov ms/it {tm.get("krylov",0)*1e3/max(cnt.get("kry_it",1),1):.3f}')
