"""minimal tab loader for run_radwave.py (replacement of the gone wt-littests vis/python common.load_tab)."""
import numpy as np, re
class _T:
    pass
def load_tab(path):
    hdr = []; t = None
    with open(path) as f:
        for ln in f:
            if not ln.startswith('#'): break
            hdr.append(ln)
    m = re.search(r'time=\s*([0-9eE.+-]+)', hdr[0]); t = float(m.group(1))
    names = hdr[-1].replace('#', '').split()
    vars_ = [n for n in names if n not in ('gid', 'i', 'j', 'k', 'x1v', 'x2v', 'x3v')]
    a = np.loadtxt(path, comments='#', ndmin=2)
    o = _T(); o.time = t; o.x1 = a[:, 2]
    off = a.shape[1] - len(vars_)
    cols = {v: a[:, off + q] for q, v in enumerate(vars_)}
    o.var = lambda n: cols[n]
    return o
