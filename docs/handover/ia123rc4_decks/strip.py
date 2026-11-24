# strip.py SUMMARY DECK... (MERGE123, COORD 11-24): drop from each deck the keys the rc3 binary never read in its own 1-cycle run
#   (orion/unread/summary.txt, job10/11: the "never read by the code" WARNING after cycle 1). Output blocks are skipped (that run
#   disabled the outputs, so their keys read as unread there). Map: deck -> run name.
import sys, re, os
runs = {}
for l in open(sys.argv[1]):
    m = re.match(r'== (\S+) rc=\S+ unread: (.*)', l)
    if m: runs[m.group(1)] = {k for k in m.group(2).split() if not k.startswith('output')}
MAP = {'IA123_t9l62': 't9l62', 'IA123_m35mms': 'm35mms', 'IA123_bsg': 'bsg', 'IA123_m35zams_wallsym': 'm35zams', 'IA123_gate_bsg_smoke': 'gate',
       'star_t9l62': 'star_t9l62', 'star_m35mms': 'star_m35mms', 'star_bsg': 'bsg', 'star3d_sp': 'star3d_sp', 'star3d_cs': 'star3d_cs'}
for d in sys.argv[2:]:
    un = runs[MAP[os.path.basename(d).split('.')[0]]]; blk = ''; out = []; n = 0
    for l in open(d):
        m = re.match(r'\s*<(\w+)>', l)
        if m: blk = m.group(1)
        k = re.match(r'\s*([A-Za-z][A-Za-z0-9_]*)\s*=', l)
        if k and not m and f'{blk}/{k.group(1)}' in un: n += 1; continue
        out.append(l)
    open(d, 'w').writelines(out); print(d, 'stripped', n)
