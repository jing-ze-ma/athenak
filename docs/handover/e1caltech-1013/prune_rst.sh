#!/bin/bash
# prodswitch-1013 (PS) prune_rst.sh RUN_DIR [KEEP=3] [RC] [JOBID]: keep only the newest KEEP restart files of a production run.
# User rule (2026-10-13): after a link that ENDED SUCCESSFULLY, if its newest rst passes a sanity check, delete the rst files older than
# the newest KEEP in RUN_DIR/rst. Never delete after a failed or STOPped link. Every action is logged to RUN_DIR/rst_prune.log.
#   success = RC == 0 (the srun rc of the link), no RUN_DIR/STOP, and this job's run log has no FATAL and no NaN in its last 60 lines.
#   sanity of the newest rst (header parsed by python, as rst_info.py does):
#     size > 0; the LEADING part (size minus the vet_rotate_avg history block, M1ROTAV1 header, which legitimately varies:
#     K-1 vs K slots, trimmed radial range) within 1 % of the next-older rst's unless that one is the cycle-0 file; a history
#     block must be consistent with its own header; time finite, > 0, and > the next-older rst's time;
#     ncycle > 0; for the rst output block of the header: 0 <= time - last_time < dt (the stored schedule is consistent).
#   order = (ncycle, name) from the headers, not the file names.
# Called from the EXIT trap that run.cfg installs (queued jobs keep the batch script they were submitted with, but they all source run.cfg).
D=${1:?}; KEEP=${2:-3}; RC=${3:-x}; JOB=${4:-${SLURM_JOB_ID:-}}
L="$D/rst_prune.log"; R="$D/rst"
log() { echo "$(date +%F_%T) job ${JOB:-?}: $*" >> "$L"; }
[ -d "$R" ] || exit 0
exec 9> "$R/.prune.lock"; flock -n 9 || { log "another prune holds the lock: skip"; exit 0; }
[ "$RC" = 0 ] || { log "link rc=$RC: no pruning"; exit 0; }
[ -f "$D/STOP" ] && { log "STOP present ($(head -c 200 "$D/STOP")): no pruning"; exit 0; }
LOGF=""; for f in "$D/run.log.$JOB" "$D/run.log"; do [ -f "$f" ] && { LOGF=$f; break; }; done
if [ -n "$LOGF" ]; then
  grep -q FATAL "$LOGF" && { log "FATAL in $LOGF: no pruning"; exit 0; }
  tail -n 60 "$LOGF" | grep -qiw nan && { log "NaN in $LOGF: no pruning"; exit 0; }
fi
python3 - "$R" "$KEEP" "$L" "${JOB:-?}" <<'PY'
import os, re, sys, glob, struct, math, time
R, keep, L, job = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
def log(s):
    open(L, 'a').write('%s job %s: %s\n' % (time.strftime('%Y-%m-%d_%H:%M:%S'), job, s))
def head(f):
    with open(f, 'rb') as fh: buf = fh.read(8000000)
    p = buf.find(b'<par_end>\n'); assert p > 0, 'no <par_end>'
    txt = buf[:p].decode('latin-1')
    o0 = p + len('<par_end>\n')
    o = o0 + 4 + 4 + 9*8 + 19*4 + 19*4
    t, dt = struct.unpack('<dd', buf[o:o+16]); nc = struct.unpack('<i', buf[o+16:o+20])[0]
    B, b = {}, None
    for ln in txt.splitlines():
        ln = ln.split('#')[0].strip(); m = re.match(r'<(\w+)>', ln)
        if m: b = m.group(1); B.setdefault(b, {}); continue
        if b and '=' in ln:
            k, v = [x.strip() for x in ln.split('=', 1)]; B[b][k] = v
    # vet_rotate_avg history (rotavg-rst-1012, marker M1ROTAV1 + uint64 nb + int64 h[nb/8]; h[1] slabs of
    # nmb_total x (nx1+2ng)(nx2+2ng)(nx3+2ng) doubles): its size legitimately varies between files (K-1 vs K slots,
    # trimmed radial range), so the size check uses only the LEADING part = file size - history bytes
    sz = os.path.getsize(f)
    nmb = struct.unpack('<i', buf[o0:o0+4])[0]
    ng, nx1, nx2, nx3 = struct.unpack('<4i', buf[o0+8+72+76:o0+8+72+76+16])
    ncell = (nx1 + 2*ng) * ((nx2 + 2*ng) if nx2 > 1 else 1) * ((nx3 + 2*ng) if nx3 > 1 else 1)
    hist, rah = 0, None
    m = buf.find(b'M1ROTAV1', o0)
    if m >= 0:
        nb = struct.unpack('<Q', buf[m+8:m+16])[0]
        assert nb % 8 == 0 and 10 <= nb // 8 <= 64, 'vet_rotate_avg header length %d' % nb
        rah = struct.unpack('<%dq' % (nb // 8), buf[m+16:m+16+nb])
        K, nsl, slots = rah[2], rah[1], rah[7]
        assert rah[0] == 1 and nsl > 0 and 2 <= K <= 64 and slots in (K - 1, K), 'vet_rotate_avg header %r' % (rah,)
        if len(rah) > 10: assert 0 <= rah[10] < nx1 + 2*ng, 'vet_rotate_avg i0 %r' % (rah[10],)
        hist = nsl * nmb * ncell * 8
        assert 0 < hist < sz, 'vet_rotate_avg history %d bytes vs file %d' % (hist, sz)
    return t, nc, B, sz - hist, hist, rah

files = sorted(glob.glob(os.path.join(R, '*.rst')))
if len(files) <= keep:
    log('%d rst file(s) <= %d: nothing to prune' % (len(files), keep)); sys.exit()
info = {}
for f in files:
    try: info[f] = head(f)
    except Exception as e: info[f] = None; log('unreadable header %s (%s)' % (os.path.basename(f), e))
good = [f for f in files if info[f] is not None]
good.sort(key=lambda f: (info[f][1], f))
new = good[-1] if good else None
def bad(msg):
    log('newest rst sanity FAILED (%s): no pruning' % msg); sys.exit()
nm = max(files, key=lambda f: (os.path.getmtime(f), f))
if info[nm] is None: bad('the most recently written file %s is unreadable' % os.path.basename(nm))
if new is None or len(good) < 2: bad('fewer than 2 readable rst files')
if nm != new: bad('the most recently written file %s is not the highest-cycle one %s' % (os.path.basename(nm), os.path.basename(new)))
t, nc, B, lead, hist, rah = info[new]; sz = os.path.getsize(new)
if sz <= 0: bad('%s size 0' % new)
if not (math.isfinite(t) and t > 0 and nc > 0): bad('%s time %r ncycle %r' % (new, t, nc))
prev = good[-2]
if info[prev][0] >= t: bad('time %r not > previous %s time %r' % (t, os.path.basename(prev), info[prev][0]))
plead = info[prev][3]
if lead <= 0: bad('leading part %d' % lead)
if info[prev][1] > 0 and abs(plead - lead) > 0.01 * plead: bad('leading part (file without the vet_rotate_avg history) %d differs > 1 %% from previous %s %d' % (lead, os.path.basename(prev), plead))
ob = [d for b, d in B.items() if b.startswith('output') and d.get('file_type') == 'rst']
for d in ob:
    dt = float(d.get('dt', '-1')); lt = float(d.get('last_time', '-1'))
    if dt > 0 and lt >= 0 and not (-1e-9 * max(1.0, t) <= t - lt < dt * (1 + 1e-9)):
        bad('rst output last_time %r inconsistent with time %r, dt %r' % (lt, t, dt))
kept = good[-keep:]
drop = [f for f in good if f not in kept]   # unreadable files are never deleted (logged above)
log('newest %s OK (t=%.6g ncycle=%d size=%d = lead %d + history %d); keep %s' % (os.path.basename(new), t, nc, sz, lead, hist, ' '.join(os.path.basename(f) for f in kept)))
for f in drop:
    s = os.path.getsize(f); os.remove(f)
    log('deleted %s (%d bytes, t=%s ncycle=%s)' % (os.path.basename(f), s, '%.6g' % info[f][0] if info[f] else '?', info[f][1] if info[f] else '?'))
PY
