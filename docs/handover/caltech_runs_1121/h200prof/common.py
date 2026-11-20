# common.py: shared helpers for picktop/gaps (orion H200PROF): kernel key as kname.py, cycle window from NewTimeStep clusters
import re, sqlite3
def open_db(p):
    db = sqlite3.connect(p); c = db.cursor(); S = dict(c.execute("select id,value from StringIds")); return c, S
def key(s):
    m = re.findall(r'(([A-Za-z_][\w:]*)\([^()]*\)::\[lambda\([^\]]*\) \(instance (\d+)\)\])', s)
    m = [x for x in m if not x[1].endswith(('par_for', 'CleanFor', 'par_reduce', 'CleanReduce')) and 'Kokkos' not in x[1]]
    if m: return m[0][1].replace('radm1::RadiationM1::', 'M1::') + '#' + m[0][2], m[0][0]
    m = re.search(r'ParallelFor<([\w:]+)', s) or re.search(r'Parallel\w+<([\w:]+)', s)
    return (m.group(1) if m else s[:60]), None
def window(c, S, c0, nc, sub):
    mk = sorted(st for st, dn in c.execute("select start,demangledName from CUPTI_ACTIVITY_KIND_KERNEL") if 'NewTimeStep' in S.get(dn, ''))
    cl = []
    for t in mk:
        if not cl or t - cl[-1][-1] > 2e8: cl.append([t])
        else: cl[-1].append(t)
    ends = [x[0] for x in cl]; j1 = min(c0 + nc, len(ends) - 1)
    return ends[c0], ends[j1], (j1 - c0)/sub


