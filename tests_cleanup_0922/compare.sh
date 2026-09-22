#!/bin/bash -l
# compare.sh <tagA> <tagB> : bitwise comparison of two gates.sh run trees.
#   .hst / .txt : compared byte for byte
#   .bin/.cbin  : compared AFTER the embedded "#--- PAR_DUMP ---" text header, whose
#                 length legitimately changes when a parameter is deleted.  The header
#                 length is the "header offset=" field of each file itself.
#   rst/        : skipped (embeds the effective parameter list)
G=/viper/u2/jinma/ATHENAK/athenak/tests_cleanup_0922
A=$G/$1; B=$G/$2
python3 - "$A" "$B" <<'PY'
import os, sys, re
a, b = sys.argv[1], sys.argv[2]
def off(p):
    with open(p, 'rb') as f:
        head = f.read(4096)
    m = re.search(rb'header offset=(\d+)', head)
    return int(m.group(1)) if m else 0
nsame = ndiff = nmiss = 0
for root, dirs, files in os.walk(a):
    dirs[:] = [d for d in dirs if d != 'rst']
    for fn in sorted(files):
        if fn == 'run.log' or fn.endswith('.log'):
            continue
        pa = os.path.join(root, fn)
        pb = os.path.join(b, os.path.relpath(pa, a))
        if not os.path.exists(pb):
            print('MISSING', os.path.relpath(pa, a)); nmiss += 1; continue
        if fn.endswith('.bin') or fn.endswith('.cbin'):
            oa, ob = off(pa), off(pb)
        else:
            oa = ob = 0
        with open(pa, 'rb') as f:
            f.seek(oa); da = f.read()
        with open(pb, 'rb') as f:
            f.seek(ob); db = f.read()
        if da == db:
            nsame += 1
        else:
            print('DIFFER', os.path.relpath(pa, a), len(da), len(db)); ndiff += 1
print('COMPARE %s vs %s : same=%d differ=%d missing=%d' % (a, b, nsame, ndiff, nmiss))
PY
