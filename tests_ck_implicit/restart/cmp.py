"""cmp.py <mode> A B.
mode rst: the restart continuation A (straight) vs B (restarted): the last rst file of
  each, bytes after <par_end>; the .hst rows of B against A's rows at the same time
  (string compare, %.17e); the bin dumps B wrote against A's of the same name.
mode run: two runs of the same thing (base vs new binary): every .hst row and every bin
  dump after its parameter text; rst files only when asked (third argument 'rst')."""
import glob
import os
import sys


def tail(f):
    d = open(f, 'rb').read()
    for m in (b'<par_end>', b'PAR_END'):
        if d.find(m) >= 0:
            return d[d.find(m):]
    return d


def hst(f):
    rows = {}
    for ln in open(f):
        if ln.startswith('#') or not ln.strip():
            continue
        rows[ln.split()[0]] = ln
    return rows


mode, a, b = sys.argv[1], sys.argv[2], sys.argv[3]
withrst = len(sys.argv) > 4 and sys.argv[4] == 'rst'
out = []
bad = 0
ha, hb = glob.glob(a + '/*.hst'), glob.glob(b + '/*.hst')
if ha and hb:
    ra, rb = hst(ha[0]), hst(hb[0])
    keys = [k for k in rb if k in ra]
    nd = sum(1 for k in keys if ra[k] != rb[k])
    if mode == 'run' and len(ra) != len(rb):
        nd += 1
    out.append('hst %d rows %s' % (len(keys), 'bitwise' if nd == 0 else '%d DIFFER' % nd))
    bad += nd
nb = ndb = 0
# bin dumps matched by content of the whole record after the parameter text (which
# carries time and cycle): the restarted run numbers its files differently
ta = set(tail(f) for f in glob.glob(a + '/bin/*.bin'))
for fb in sorted(glob.glob(b + '/bin/*.bin')):
    tb = tail(fb)
    nb += 1
    if tb not in ta:
        ndb += 1
out.append('bin %d %s' % (nb, 'bitwise' if ndb == 0 else '%d DIFFER' % ndb))
bad += ndb
if mode == 'rst' or withrst:
    sa = sorted(glob.glob(a + '/rst/*.rst'))
    sb = sorted(glob.glob(b + '/rst/*.rst'))
    if sa and sb:
        same = tail(sa[-1]) == tail(sb[-1])
        out.append('rst %s/%s %s' % (os.path.basename(sa[-1]), os.path.basename(sb[-1]),
                                     'bitwise' if same else 'DIFFER'))
        bad += 0 if same else 1
    else:
        out.append('rst MISSING')
        bad += 1
print('%-60s %s  => %s' % (b.split('ckrst_0924/')[-1], '; '.join(out),
                           'BITWISE' if bad == 0 else 'NOT BITWISE'))
