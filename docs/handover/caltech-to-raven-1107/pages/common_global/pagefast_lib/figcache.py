"""pagefast-1102: exact caches for the RSG/HG pages' mk_figs.py (both import this module).
lazy_npz(files, cache_dir): list of read-only mappings that serve the 0-d members (t, Rph_med, Teff_med, nonfinite, vet, ver, ...)
  and .files from a per-dump cache (keyed by name, size, mtime_ns), opening the real npz only when an array is asked for.
  The cached values are the npz's own 0-d arrays (pickled), so float()/int()/bool() of them are bit-identical.
logscan(texts_files): per-file cached regex results for the status block (keyed by path, size, mtime_ns), so only a new or
  growing log is scanned again."""
import os, re, pickle, zipfile
from collections.abc import Mapping
import numpy as np


class _Lazy(Mapping):
    def __init__(self, f, ent):
        self._f, self._e, self._z = f, ent, None

    @property
    def files(self):
        return list(self._e['files'])

    def _real(self):
        if self._z is None:
            self._z = np.load(self._f)
        return self._z

    def __getitem__(self, k):
        s = self._e['scal']
        if k in s:
            return s[k].copy()
        if k not in self._e['files']:
            raise KeyError(k)
        return self._real()[k]

    def __iter__(self):
        return iter(self._e['files'])

    def __len__(self):
        return len(self._e['files'])

    def __contains__(self, k):
        return k in self._e['files']


def _entry(f):
    z = np.load(f)
    zf = zipfile.ZipFile(f)
    small = {i.filename[:-4] for i in zf.infolist() if i.file_size < 4096}
    zf.close()
    scal = {}
    for k in z.files:
        if k in small:
            a = z[k]
            if isinstance(a, np.ndarray) and a.ndim == 0:
                scal[k] = a
    return dict(files=list(z.files), scal=scal)


def lazy_npz(files, cache_dir):
    cf = os.path.join(cache_dir, '.figcache_scal.pkl')
    try:
        C = pickle.load(open(cf, 'rb'))
    except Exception:
        C = {}
    out, new, chg = [], {}, False
    for f in files:
        st = os.stat(f); k = os.path.basename(f); sig = (st.st_size, st.st_mtime_ns)
        e = C.get(k)
        if e is None or e['sig'] != sig:
            e = _entry(f); e['sig'] = sig; chg = True
        new[k] = e
        out.append(_Lazy(f, e))
    if chg or len(new) != len(C):
        try:
            tmp = cf + '.%d.tmp' % os.getpid()
            pickle.dump(new, open(tmp, 'wb'), protocol=4); os.replace(tmp, cf)
        except Exception:
            pass
    return out


PATS = {'el': re.compile(r'elapsed=([0-9.eE+-]+) cycle=([0-9]+) time=([0-9.eE+-]+)'), 'FATAL': re.compile(r'FATAL'),
        'DIVERGED': re.compile(r'DIVERGED'), 'nan': re.compile(r'\bnan\b|\bNaN\b'), 'nc': re.compile(r'NON-CONVERGED=[1-9]')}
LCACHE = os.environ.get('LOGCACHE', '/resnick/groups/carnegie_poc/jingze/pagefast_lib/logcache.pkl')


def logscan(logs, lks, rd):
    """-> el (the last <= 20 elapsed matches of every log, in order: el[-1], el[-min(len, 20)] and len(el) > 3 are
    the same as for the full list), counts {FATAL, DIVERGED, nan, nc} over logs + lks, as the old findall on the joined text."""
    try:
        C = pickle.load(open(LCACHE, 'rb'))
    except Exception:
        C = {}
    chg = False
    el, cnt = [], dict(FATAL=0, DIVERGED=0, nan=0, nc=0)
    for f, is_log in [(x, True) for x in logs] + [(x, False) for x in lks]:
        try:
            st = os.stat(f); sig = (st.st_size, st.st_mtime_ns)
        except OSError:
            sig = None
        key = (f, is_log)
        e = C.get(key)
        if sig is None or e is None or e['sig'] != sig:
            txt = rd(f)
            e = dict(sig=sig, el=PATS['el'].findall(txt)[-20:] if is_log else [],
                     **{k: len(PATS[k].findall(txt)) for k in cnt})
            if sig is not None:
                C[key] = e; chg = True
        el += e['el']
        for k in cnt:
            cnt[k] += e[k]
    if chg:
        try:
            tmp = LCACHE + '.%d.tmp' % os.getpid()
            pickle.dump(C, open(tmp, 'wb'), protocol=4); os.replace(tmp, LCACHE)
        except Exception:
            pass
    return el, cnt
