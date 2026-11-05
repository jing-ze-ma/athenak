"""binprec.py (USER 11-05: the float32-dump trap is structural): gate comparisons must use DOUBLE output.
require_double(paths, allow=False, what='gate compare'): reads each AthenaK .bin header ('size of variable=N') and
exits(3) with a clear message if any file is float32, unless allow (--allow-float32) is set.
require_tab_digits(paths, allow=False): the same for .tab dumps (needs >= 16 significant digits per value).
Copies: tools/binprec.py (scripts outside the repo) and tst/layout_invariance/binprec.py (repo)."""
import sys, re

def bin_var_size(path):
    with open(path, 'rb') as f:
        head = f.read(4096).decode(errors='replace')
    m = re.search(r'size of variable\s*=\s*(\d+)', head)
    return int(m.group(1)) if m else None

def _die(msg):
    sys.stderr.write('### REFUSED: ' + msg + '\n'); print('REFUSED ' + msg); sys.exit(3)

def require_double(paths, allow=False, what='gate compare'):
    if allow: return
    bad = [p for p in paths if bin_var_size(p) != 8]
    if bad:
        _die('%s on float32 (or unreadable) .bin dumps: %s%s -- write them with <outputN>/bin_double = true '
             '(float32 rounding hides/creates differences at 1e-8..1e-7 relative), or pass --allow-float32'
             % (what, bad[0], ' (+%d more)' % (len(bad) - 1) if len(bad) > 1 else ''))

def require_tab_digits(paths, allow=False, what='tab compare'):
    if allow: return
    for p in paths:
        for l in open(p):
            if l.startswith('#') or not l.strip(): continue
            v = l.split()[-1]; mant = re.sub(r'[eE].*$', '', v).lstrip('+-').replace('.', '').lstrip('0')
            if len(mant) < 16:
                _die('%s on %s: values carry %d significant digits (need >= 16: data_format = %%.17e), or pass --allow-float32'
                     % (what, p, len(mant)))
            break

def pop_allow(argv):
    """remove --allow-float32 from argv (in place); return True if it was there"""
    if '--allow-float32' in argv:
        argv.remove('--allow-float32'); return True
    return False
