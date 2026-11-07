#!/usr/bin/env python3
"""preflight.py [--md5 MD5] BIN athena-args...   (run on a LOGIN node before submitting / before a job can start)
Mirrors ParameterInput::ModifyFromCmdline (prod-1107v): fresh (-i) -> every block/key on the command line must exist in
the deck; restart (-r) -> the block must exist in the rst's input (absent keys are added). Also: deck/rst/BIN exist, BIN md5,
every *_file / *_table value (deck or cmdline) is a readable file, -t is HH:MM:SS with no negatives. Exit 1 on any problem."""
import sys, os, re, hashlib
a = sys.argv[1:]; md5 = None
if a and a[0] == '--md5': md5, a = a[1], a[2:]
binp, args = a[0], a[1:]; errs = []
def blocks_of(text):
    B, cur = {}, None
    for l in text.split('\n'):
        s = l.split('#', 1)[0].strip()
        m = re.match(r'<\s*([^>]+?)\s*>', s)
        if m: cur = m.group(1); B.setdefault(cur, {}); continue
        if cur and '=' in s: k, v = s.split('=', 1); B[cur][k.strip()] = v.strip()
    return B
if not os.access(binp, os.X_OK): errs.append('binary missing/not executable: ' + binp)
elif md5:
    h = hashlib.md5(open(binp, 'rb').read()).hexdigest()
    if h != md5: errs.append('binary md5 %s != %s' % (h, md5))
mode = path = None
for i, x in enumerate(args):
    if x in ('-i', '-r') and i + 1 < len(args): mode, path = x, args[i + 1]
    if x == '-t' and i + 1 < len(args) and not re.fullmatch(r'\d+:\d{2}:\d{2}', args[i + 1]): errs.append('bad -t ' + args[i + 1])
if not path or not os.path.isfile(path): errs.append('input %s missing: %s' % (mode, path))
else:
    if mode == '-i': text = open(path, errors='replace').read()
    else:
        raw = open(path, 'rb').read(4 << 20); e = raw.find(b'<par_end>')
        if e < 0: errs.append('no <par_end> in first 4 MB of rst'); text = ''
        else: text = raw[:e].decode(errors='replace')
    B = blocks_of(text)
    for x in args:
        m = re.fullmatch(r'([^/=\s]+)/([^=\s]+)=(.*)', x)
        if not m: continue
        b, k, v = m.groups()
        if b not in B: errs.append('block <%s> not in %s (%s)' % (b, mode, x)); continue
        if mode == '-i' and k not in B[b]: errs.append('%s/%s not in deck: fresh start will FATAL' % (b, k))
        B[b][k] = v
    for b in B:
        for k, v in B[b].items():
            if re.search(r'(_file|_table)$', k) and v.startswith('/') and not os.path.isfile(v):
                errs.append('%s/%s file missing: %s' % (b, k, v))
for e in errs: print('PREFLIGHT FAIL:', e)
print('PREFLIGHT', 'OK' if not errs else 'FAILED (%d)' % len(errs), '|', mode, path)
sys.exit(1 if errs else 0)
