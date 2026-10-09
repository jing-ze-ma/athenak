#!/usr/bin/env python3
# rst_info.py RST [OVERLAY]: prints "time ncycle tlim" then a line of outputN/last_time keys
# = floor(time/dt)*dt for every <outputN> with dt > 0 (rst header, overlay adds/overrides).
import re, struct, sys, math
f = open(sys.argv[1], 'rb'); buf = f.read(2000000); f.close()
p = buf.find(b'<par_end>\n'); assert p > 0
txt = buf[:p].decode('latin-1')
o = p + len('<par_end>\n') + 4 + 4 + 9*8 + 19*4 + 19*4
t, dt = struct.unpack('<dd', buf[o:o+16]); nc = struct.unpack('<i', buf[o+16:o+20])[0]
def blocks(s):
    b, d = None, {}
    for ln in s.splitlines():
        ln = ln.split('#')[0].strip()
        m = re.match(r'<(\w+)>', ln)
        if m: b = m.group(1); d.setdefault(b, {}); continue
        if b and '=' in ln:
            k, v = [x.strip() for x in ln.split('=', 1)]; d[b][k] = v
    return d
B = blocks(txt)
if len(sys.argv) > 2:
    for k, v in blocks(open(sys.argv[2]).read()).items(): B.setdefault(k, {}).update(v)
keys = []
for b, d in sorted(B.items()):
    if b.startswith('output') and float(d.get('dt', '-1')) > 0:
        odt = float(d['dt']); keys.append('%s/last_time=%.17g' % (b, math.floor(t/odt)*odt))
print('%.17g %d %s' % (t, nc, B['time']['tlim']))
print(' '.join(keys))
