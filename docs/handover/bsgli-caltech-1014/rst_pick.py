#!/usr/bin/env python3
# prodswitch-1013 (PS) rst_pick.py RSTDIR [--quarantine]: pick the newest VALID restart file of a production run dir.
# Output (stdout) when one exists, else nothing (no rst at all); exit 2 when rst files exist but none is valid:
#   line 1: "<path> <time> <ncycle> <tlim>"
#   line 2: outputN/last_time = floor(time/dt)*dt for every <outputN> with dt > 0 (as rst_info.py: the output cadence is kept)
# Valid = the header parses (as rst_info.py / prune_rst.sh do), time finite and >= 0, a vet_rotate_avg history block (M1ROTAV1)
#   consistent with its own header and inside the file, and, if an older readable rst with ncycle > 0 exists, the LEADING part
#   (size minus the history, which legitimately varies) within 1 % of that file's (catches a truncated / half-written file).
# Order: (ncycle, name) from the headers. Files newer than the pick that fail are reported on stderr; with --quarantine they are
#   moved to RSTDIR/invalid/ (so a link that takes the newest *.rst by name restarts from the pick). Nothing is ever deleted.
import os, re, sys, glob, struct, math
R = sys.argv[1]; quar = '--quarantine' in sys.argv[2:]
def head(f):
    with open(f, 'rb') as fh: buf = fh.read(8000000)
    p = buf.find(b'<par_end>\n'); assert p > 0, 'no <par_end>'
    txt = buf[:p].decode('latin-1'); o0 = p + len('<par_end>\n')
    o = o0 + 4 + 4 + 9*8 + 19*4 + 19*4
    assert len(buf) >= o + 20, 'header truncated'
    t, dt = struct.unpack('<dd', buf[o:o+16]); nc = struct.unpack('<i', buf[o+16:o+20])[0]
    B, b = {}, None
    for ln in txt.splitlines():
        ln = ln.split('#')[0].strip(); m = re.match(r'<(\w+)>', ln)
        if m: b = m.group(1); B.setdefault(b, {}); continue
        if b and '=' in ln:
            k, v = [x.strip() for x in ln.split('=', 1)]; B[b][k] = v
    sz = os.path.getsize(f)
    nmb = struct.unpack('<i', buf[o0:o0+4])[0]
    ng, nx1, nx2, nx3 = struct.unpack('<4i', buf[o0+156:o0+172])
    ncell = (nx1 + 2*ng) * ((nx2 + 2*ng) if nx2 > 1 else 1) * ((nx3 + 2*ng) if nx3 > 1 else 1)
    hist = 0; m = buf.find(b'M1ROTAV1', o0)
    if m >= 0:
        nb = struct.unpack('<Q', buf[m+8:m+16])[0]
        assert nb % 8 == 0 and 10 <= nb // 8 <= 64, 'vet_rotate_avg header length %d' % nb
        h = struct.unpack('<%dq' % (nb // 8), buf[m+16:m+16+nb])
        K, nsl, slots = h[2], h[1], h[7]
        assert h[0] == 1 and nsl > 0 and 2 <= K <= 64 and slots in (K - 1, K), 'vet_rotate_avg header %r' % (h,)
        if len(h) > 10: assert 0 <= h[10] < nx1 + 2*ng, 'vet_rotate_avg i0 %r' % (h[10],)
        hist = nsl * nmb * ncell * 8
        assert 0 < hist < sz - (m + 16 + nb), 'vet_rotate_avg history %d bytes vs file %d' % (hist, sz)
    assert math.isfinite(t) and t >= 0 and nc >= 0 and math.isfinite(dt) and dt > 0, 'time %r dt %r ncycle %r' % (t, dt, nc)
    return dict(t=t, nc=nc, B=B, lead=sz - hist)
files = sorted(glob.glob(os.path.join(R, '*.rst')))
info = {}
for f in files:
    try: info[f] = head(f)
    except Exception as e: info[f] = None; print('rst_pick: unreadable %s (%s)' % (f, e), file=sys.stderr)
good = sorted([f for f in files if info[f] is not None], key=lambda f: (info[f]['nc'], f))
pick = None
for i in range(len(good) - 1, -1, -1):
    f = good[i]; ref = [g for g in good[:i] if info[g]['nc'] > 0]
    if ref and abs(info[f]['lead'] - info[ref[-1]]['lead']) > 0.01 * info[ref[-1]]['lead']:
        print('rst_pick: %s leading part %d vs %s %d (> 1 %%): rejected' % (f, info[f]['lead'], ref[-1], info[ref[-1]]['lead']), file=sys.stderr)
        continue
    pick = f; break
if files and pick is None:   # rst files exist but none is valid: never start fresh over them
    print('rst_pick: %d rst file(s), none valid' % len(files), file=sys.stderr); sys.exit(2)
bad = [f for f in files if f > pick] if pick else []   # newer by name than the pick: unreadable or rejected
for f in bad: print('rst_pick: %s is newer than the pick %s but invalid' % (f, pick), file=sys.stderr)
if quar and bad:
    os.makedirs(os.path.join(R, 'invalid'), exist_ok=True)
    for f in bad:
        os.rename(f, os.path.join(R, 'invalid', os.path.basename(f)))
        print('rst_pick: quarantined %s -> invalid/' % f, file=sys.stderr)
if pick is not None:
    d = info[pick]; B = d['B']; t = d['t']; keys = []
    for b, v in sorted(B.items()):
        if b.startswith('output') and float(v.get('dt', '-1')) > 0:
            odt = float(v['dt']); keys.append('%s/last_time=%.17g' % (b, math.floor(t/odt)*odt))
    print('%s %.17g %d %s' % (pick, t, d['nc'], B['time']['tlim']))
    print(' '.join(keys))
