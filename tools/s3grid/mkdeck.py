"""mkdeck.py (s3grid-1017): child decks from the S3 parent parameter dump (rst 00046) with the extended fine radial grid.
usage: mkdeck.py PARENT_PAR FIT_TXT WB_RMAX_RPH OUTDIR
writes OUTDIR/s3g_tmpl.athinput (fresh template run: he_ic_equil on the new grid, 1 cycle, rst every cycle) and
OUTDIR/s3g_run.athinput (the run deck; remap_r.py adds start_time / he_ic_remap_file / output schedules)."""
import re
import sys

RPH = 2.7025e13
par, fit, wbr, od = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
txt = open(par).read()
txt = txt[:txt.find('<par_end>')]
kv = {}
for ln in open(fit):
    if ln.startswith('#') or '=' not in ln:
        continue
    k, v = [s.strip() for s in ln.split('=')]
    kv[k] = v
grid = {'nx1': kv['nx1']}
for k in kv:
    if k.startswith('f_stretch_r_'):
        grid[k] = kv[k]


def setkey(t, blk, key, val, cmt):
    lines = t.split('\n')
    bi = None
    for n, ln in enumerate(lines):
        if re.match(r'\s*<\s*%s\s*>' % re.escape(blk), ln):
            bi = n
            continue
        if bi is not None and re.match(r'\s*<', ln):
            break
        if bi is not None and re.match(r'\s*%s\s*=' % re.escape(key), ln):
            lines[n] = '%-26s = %s   # %s' % (key, val, cmt)
            return '\n'.join(lines)
    if bi is None:
        raise SystemExit('no block ' + blk)
    lines.insert(bi + 1, '%-26s = %s   # %s' % (key, val, cmt))
    return '\n'.join(lines)


def delkey(t, blk, key):
    lines = t.split('\n')
    out, inb = [], False
    for ln in lines:
        m = re.match(r'\s*<\s*([^>]+)\s*>', ln)
        if m:
            inb = (m.group(1).strip() == blk)
        if inb and re.match(r'\s*%s\s*=' % re.escape(key), ln):
            continue
        out.append(ln)
    return '\n'.join(out)


C = 's3grid-1017 (gridfit2.py: S3 faces to 1.08 R_ph, dr/r 0.3 % to 1.72 R_ph, then x1.04/cell)'
t = txt
for k, v in grid.items():
    t = setkey(t, 'mesh', k, v, C)
t = setkey(t, 'meshblock', 'nx1', grid['nx1'], C)
t = setkey(t, 'hydro', 'wb_rmax', '%.6e' % (wbr * RPH), 's3grid-1017: x1 WB to %.2f R_ph (was 1.1)' % wbr)
t = t.replace('<comment>\n', '<comment>\ns3grid-1017: S3 (AG Car A 268x32x32 test) with the fine radial zone and wb_rmax extended to '
              '%.2f R_ph (nx1 %s); BGDIAG grid fix test. ' % (wbr, grid['nx1']), 1)
for blk in ('output1', 'output2', 'output3', 'output4', 'output5'):
    t = delkey(t, blk, 'last_time')
    t = delkey(t, blk, 'file_number')
# run deck: long-stat outputs (as rol2u_safety_1015 vlong: hst 1e4, hydro_w 2e6, m1 off, log 1e5, rst 1e7), ndiag 1
r = t
r = setkey(r, 'time', 'tlim', '3.3e8', 's3grid-1017 run (>= 2.5e8 + margin; restarted links)')
r = setkey(r, 'time', 'ndiag', '1', 's3grid-1017')
r = setkey(r, 'output2', 'dt', '2e6', 's3grid-1017 (as vlong)')
r = setkey(r, 'output3', 'dt', '1e30', 's3grid-1017 (m1 dumps off, as vlong)')
r = setkey(r, 'output5', 'dt', '1e7', 's3grid-1017 (as vlong)')
open(od + '/s3g_run.athinput', 'w').write(r + '<par_end>\n' if False else r)
# template: fresh, 1 cycle, rst after it, nothing else
m = t
m = setkey(m, 'time', 'tlim', '3.0e7', 's3grid-1017 template')
m = setkey(m, 'time', 'nlim', '1', 's3grid-1017 template: 1 cycle')
m = setkey(m, 'output2', 'dt', '1e30', 'template')
m = setkey(m, 'output3', 'dt', '1e30', 'template')
m = setkey(m, 'output5', 'dt', '1e30', 'template')
m = setkey(m, 'output5', 'dcycle', '1', 'template: rst after cycle 1')
open(od + '/s3g_tmpl.athinput', 'w').write(m)
print('wrote', od + '/s3g_run.athinput', od + '/s3g_tmpl.athinput', 'nx1', grid['nx1'])
