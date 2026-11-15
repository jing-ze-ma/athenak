"""ana_pos.py: breakdown of pos.out violations per variant by test, cfl, magnitude."""
import re, collections
cur, data = None, collections.OrderedDict()
for ln in open(__import__('sys').argv[1] if len(__import__('sys').argv) > 1 else 'pos.out'):
    m = re.match(r'(\S+)\s+cases', ln)
    if m:
        cur = m.group(1); data[cur] = []; continue
    m = re.match(r'\s+(\w+)\s+chi dx (\S+)\s+(\w+)\s+cfl (\S+)\s+minE (\S+) max\(\|F\|-cE\) (\S+)', ln)
    if m and cur:
        t, x, k, c, e, f = m.groups(); data[cur].append((t, float(x), k, float(c), -float(e), float(f)))
for v, rows in data.items():
    w = [max(r[4], r[5]) for r in rows]
    bt = collections.Counter(r[0] for r in rows)
    bc = collections.Counter(r[3] for r in rows)
    big = sum(1 for x in w if x > 1e-3); mid = sum(1 for x in w if 1e-6 < x <= 1e-3)
    print('%-30s n %3d  >1e-3 %3d  1e-6..1e-3 %3d  <1e-6 %3d | %s | cfl %s' % (v, len(rows), big, mid,
          len(rows) - big - mid, ' '.join('%s:%d' % kv for kv in sorted(bt.items())),
          ' '.join('%g:%d' % kv for kv in sorted(bc.items()))))
