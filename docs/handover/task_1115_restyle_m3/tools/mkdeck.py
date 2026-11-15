# mkdeck.py <deck> <out> [block/par=value ...]: copy deck, apply the overrides (adding
# missing parameters/blocks), append an hst block (every cycle) unless one exists, and bin
# dumps of the evolved variables every 5 cycles.
import sys, re
s = open(sys.argv[1]).read()
lines = s.split('\n')
for kv in sys.argv[3:]:
    key, val = kv.split('=', 1)
    blk, par = key.split('/', 1)
    cur = None; done = False; bend = None
    for n, l in enumerate(lines):
        m = re.match(r'^\s*<([^>]+)>', l)
        if m:
            if cur == blk and bend is None: bend = n
            cur = m.group(1).strip(); continue
        if cur == blk and re.match(r'^\s*%s\s*=' % re.escape(par), l):
            lines[n] = '%s = %s' % (par, val); done = True
    if not done:
        idx = [n for n, l in enumerate(lines) if re.match(r'^\s*<%s>' % re.escape(blk), l)]
        if idx:
            lines.insert(idx[0] + 1, '%s = %s' % (par, val))
        else:
            lines += ['', '<%s>' % blk, '%s = %s' % (par, val)]
s = '\n'.join(lines)
blocks = set(b.strip() for b in re.findall(r'^\s*<([^>]+)>', s, re.M))
n = max([int(m) for m in re.findall(r'^\s*<output(\d+)>', s, re.M)] + [0])
types = re.findall(r'^\s*file_type\s*=\s*(\w+)', s, re.M)
out = '\n'
vs = []
if 'hydro' in blocks: vs.append('hydro_u')
if 'mhd' in blocks: vs += ['mhd_u', 'mhd_bcc']
if 'rad_m1' in blocks: vs.append('m1')
if 'hst' not in types:
    n += 1; out += '<output%d>\nfile_type = hst\ndcycle = 1\n\n' % n
for k, v in enumerate(vs):
    n += 1
    out += '<output%d>\nfile_type = bin\nvariable = %s\ndcycle = 5\nid = cmp%d\n\n' % (n, v, k)
open(sys.argv[2], 'w').write(s + out)
