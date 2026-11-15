"""idort-1019 shadow trees for the v10 page scripts (extract.py / plot.py):
cart: p19runs/cart/<case>/FIN_<x>   -> idort1019/runs/cart/<case>/FIN3_<x>/p   (PROD, ab758e94)
      p19runs/cart/<case>/PRODd_<x> -> fin3_1012/runs_281dd686/<case>/FIN3_<x> (v10 page data, FIN3 281dd686) = 'before'
sp:   p19runs/sp/<case>/<arm>_<n>   -> idort1019/runs/sp/<case>/<arm>_<hash>_<n>/p  (own runs only; p2 likewise)
lad:  p19runs/lad/<case>/<dir>      -> idort1019/runs/lad/<case>/<dir>/p
Only arms with rc.txt == 0 are linked."""
import os, glob, re, shutil
R = '/raven/ptmp/jinma/idort1020/runs'; OUT = os.environ.get('OUT19', '/raven/ptmp/jinma/idort1020/page/p19runs'); F3 = '/raven/ptmp/jinma/fin3_1012/runs_281dd686'
ARM = os.environ.get('ARM', 'p')
shutil.rmtree(OUT, ignore_errors=True); n = 0
def ok(d):
    try: return open(d + '/rc.txt').read().strip() == '0'
    except Exception: return False
def ln(src, dst):
    global n
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.lexists(dst): os.remove(dst)
    os.symlink(src, dst); n += 1
for d in sorted(glob.glob(R + '/cart/*/FIN3_*/' + ARM)):
    case, run = d.split('/')[-3:-1]
    if ok(d): ln(d, '%s/cart/%s/FIN_%s' % (OUT, case, run[5:]))
for d in []:   # idort-1020: single PROD arm, no v10 overlay
    case, run = d.split('/')[-2:]
    if ok(d): ln(os.path.realpath(d), '%s/cart/%s/PRODd_%s' % (OUT, case, run[5:]))
for g in ('sp', 'p2'):
    for d in sorted(glob.glob(R + '/%s/*/*/%s' % (g, ARM))):
        case, run = d.split('/')[-3:-1]
        m = re.match(r'^(.+)_([0-9a-f]{6})_(\d+)$', run)
        if m and ok(d): ln(d, '%s/%s/%s/%s_%s' % (OUT, g, case, m.group(1), m.group(3)))
for d in sorted(glob.glob(R + '/lad/*/*/' + ARM)):
    case, run = d.split('/')[-3:-1]
    if ok(d): ln(d, '%s/lad/%s/%s' % (OUT, case, run))
print(n, 'links')
