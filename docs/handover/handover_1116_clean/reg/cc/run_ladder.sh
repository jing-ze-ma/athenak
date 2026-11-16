#!/bin/bash
# run_ladder.sh ARM BIN "rad_m1 lines" DECK FACTOR TLIM
# One run of DECK at FACTOR x its resolution (mesh = one MeshBlock), to time TLIM,
# in <deck>/<ARM>_x<FACTOR>.
ARM=$1; B=$2; X=$3; DK=$4; F=$5; T=$6
H=$(cd "$(dirname "$0")" && pwd); SRC=$H/decks/$DK.athinput
D=$H/suite/$DK/${ARM}_x$F; mkdir -p $D; cd $D || exit 1
python3 - "$X" $SRC in.athinput $F $T <<'PY'
import re, sys
x, src, dst, f, t = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
keys = [k.split('=')[0].strip() for k in x.split(';') if '=' in k]
src_lines = [ln for ln in open(src).read().splitlines()
             if not any(re.match(r'\s*%s\s*=' % k, ln) for k in keys)]
s = '\n'.join(src_lines).replace('<rad_m1>', '<rad_m1>\n' + x.replace(';', '\n'), 1)
out, blk = [], ''
for ln in s.splitlines():
    m = re.match(r'\s*<(\w+)>', ln)
    if m: blk = m.group(1)
    m = re.match(r'(\s*)(nx1|nx2)(\s*=\s*)(\d+)(.*)', ln)
    if m and blk in ('mesh', 'meshblock') and int(m.group(4)) > 4:
        ln = '%s%s%s%d%s' % (m.group(1), m.group(2), m.group(3), int(int(m.group(4))*f), m.group(5))
    m = re.match(r'(\s*)(nlim|tlim)(\s*=\s*)(\S+)(.*)', ln)
    if m and blk == 'time':
        ln = '%s%s%s%s%s' % (m.group(1), m.group(2), m.group(3), '-1' if m.group(2) == 'nlim' else t, m.group(5))
    out.append(ln)
    if re.match(r'\s*<output\d+>', ln):
        out.append('bin_double = true')
open(dst, 'w').write('\n'.join(out) + '\n')
PY
/usr/bin/time -f "wall %e" nice -n 5 timeout 28800 $B -i in.athinput > run.log 2>&1; echo "rc=$?" >> run.log
