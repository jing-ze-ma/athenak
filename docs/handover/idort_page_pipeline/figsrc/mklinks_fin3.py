"""FIN3 final page: shadow run trees so extract.py reads the final-scheme runs under the old names.
cart: fin3runs/cart/<case>/FIN_<x>  -> fin3_1012/runs_281dd686/<case>/FIN3_<x>   (IN's final-config reruns, frozen=false)
      fin3runs/cart/<case>/PRODd_<x> -> fin3_1012/runs_281dd686/<case>/F2O_<x>   (FIN2 control = 'before')
sp:   fin3runs/sp/<case>/<arm>_<n>   -> fin3_1012/sp_FIN3_281dd686/<case>/<arm>+<FIN3 keys>_<n>
Only finished runs (rc.txt == 0 for Cartesian) are linked; re-run this as runs land."""
import os, re, glob
F3 = '/raven/ptmp/jinma/fin3_1012'
OUT = '/raven/ptmp/jinma/idort_page/fin3runs'
n = 0
for d in sorted(glob.glob(F3 + '/runs_281dd686/*/*')):
    case, run = d.split('/')[-2:]
    if run.startswith('FIN3_'):
        new = 'FIN_' + run[5:]
        try:
            if open(d + '/rc.txt').read().strip() != '0':
                continue
        except Exception:
            continue
    elif run.startswith('F2O_'):
        new = 'PRODd_' + run[4:]
    else:
        continue
    os.makedirs(OUT + '/cart/' + case, exist_ok=True)
    L = OUT + '/cart/%s/%s' % (case, new)
    if os.path.islink(L):
        os.remove(L)
    os.symlink(os.path.realpath(d), L); n += 1
for d in sorted(glob.glob(F3 + '/sp_FIN3_281dd686/*/*')):
    case, run = d.split('/')[-2:]
    m = re.match(r'^([^+]+)\+.*_(\d+)$', run)
    if not m:
        continue
    os.makedirs(OUT + '/sp/' + case, exist_ok=True)
    L = OUT + '/sp/%s/%s_%s' % (case, m.group(1), m.group(2))
    if os.path.islink(L):
        os.remove(L)
    os.symlink(os.path.realpath(d), L); n += 1
print(n, 'links')
