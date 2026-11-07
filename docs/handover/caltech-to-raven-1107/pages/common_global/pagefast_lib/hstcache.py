"""pagefast-1102: incremental, exact cache for the pages' readhst() (hst files of 100-250 MB were re-parsed on every refresh).
readhst(f) returns exactly what the original readhst returned: labels from the first '#  [' line, every non-'#' line split on
whitespace, rows with the label count kept, float64 via np.array(rows, float), then the backward 'keep strictly decreasing
time' restart filter.  The parsed complete lines are cached (CACHE/<md5 of path>.npy + .json) with the byte offset and an md5 of
the bytes up to it; a call re-parses only the bytes after the offset (the prefix md5 must still match, else a full re-parse:
handles rewrites/rewinds).  A trailing line without newline is parsed for this call but never cached."""
import os, json, hashlib
import numpy as np
CACHE = os.environ.get('HSTCACHE', '/resnick/groups/carnegie_poc/jingze/pagefast_lib/hstcache')


def _md5(fh, n):
    h = hashlib.md5(); fh.seek(0); left = n
    while left > 0:
        b = fh.read(min(left, 1 << 24))
        if not b:
            break
        h.update(b); left -= len(b)
    return h


def _parse(txt, nc):
    rows = [ln.split() for ln in txt.splitlines() if not ln.startswith('#')]
    rows = [x for x in rows if len(x) == nc]
    return np.array(rows, float).reshape(len(rows), nc)


def _lab(fh):
    fh.seek(0)
    for raw in fh:
        ln = raw.decode()
        if ln.startswith('#  ['):
            return [x.split('=')[1] for x in ln[1:].split()]
    return None


def readhst_raw(f):
    os.makedirs(CACHE, exist_ok=True)
    key = hashlib.md5(os.path.abspath(f).encode()).hexdigest()
    cj, cn = os.path.join(CACHE, key + '.json'), os.path.join(CACHE, key + '.npy')
    with open(f, 'rb') as fh:
        size = os.fstat(fh.fileno()).st_size
        meta = None
        try:
            meta = json.load(open(cj))
            h = _md5(fh, meta['off']) if meta['off'] <= size else None
            if h is None or h.hexdigest() != meta['md5']:
                meta = None
            else:
                d0 = np.load(cn)
                if d0.shape != (meta['nrow'], len(meta['lab'])):
                    meta = None
        except Exception:
            meta = None
        if meta is None:
            lab = _lab(fh)
            nc = len(lab)
            off = 0; d0 = np.zeros((0, nc)); h = hashlib.md5()
        else:
            lab = meta['lab']; nc = len(lab); off = meta['off']
        fh.seek(off)
        tail = fh.read()
    cut = tail.rfind(b'\n') + 1          # complete lines only go into the cache
    dnew = _parse(tail[:cut].decode(), nc)
    dpart = _parse(tail[cut:].decode(), nc)
    d = np.concatenate([d0, dnew]) if len(dnew) else d0
    if len(dnew) or meta is None:
        noff = off + cut
        h.update(tail[:cut]); m = h.hexdigest()
        tn = cn + '.%d.tmp.npy' % os.getpid(); tj = cj + '.%d.tmp' % os.getpid()
        np.save(tn, d); json.dump(dict(off=noff, md5=m, nrow=len(d), lab=lab), open(tj, 'w'))
        os.replace(tn, cn); os.replace(tj, cj)
    if len(dpart):
        d = np.concatenate([d, dpart])
    return lab, d


def readhst(f):
    lab, d = readhst_raw(f)
    t = d[:, 0]
    keep = np.zeros(len(t), bool)
    tmin = np.inf
    for i in range(len(t) - 1, -1, -1):
        if t[i] < tmin:
            keep[i] = True
            tmin = t[i]
    d = d[keep]
    return {n: d[:, i] for i, n in enumerate(lab)}
