"""ppcache (pagefast-1102): shared fast reader + per-dump cache for AthenaK .bin dumps (formats 1.1 / 1.2).

API
  read_binary(path, vars=None)       drop-in for vis/python bin_convert.read_binary: same dict keys and the SAME values (the
                                     dump's own float32 bytes), but only the variables in `vars` are read (np.memmap of the
                                     per-MeshBlock records; no full-file read when a subset is asked for).
  load(path, fields, cache=True)     {field: array (nmb, n3, n2, n1) float32 of the dump} + '_hdr' (time, cycle, geometry,
                                     logical, coords, var list).  With cache=True each field is also kept as a contiguous
                                     .npy under CACHE/<run tag>/<dump name>/ (read back with mmap: a later tool that asks for
                                     one field reads only that field).  Raw fields are exact copies (no precision change).
  derived(path, name, version, fn, deps, dtype=None)
                                     cached derived field: fn({dep: array}) -> array (or dict of arrays); stored once per
                                     (dump, name, version); later calls by any tool return the stored result.  Use a new
                                     version string whenever fn changes.  float64 is kept unless the caller passes dtype.
  geometry(path)                     per-run cached header/geometry (mb_logical, mb_geometry, coords, mesh) without data.
  prune(max_gb=300)                  drop the least recently used dump entries beyond max_gb (raw fields are dump-sized).
Cache root: $PPCACHE or /resnick/groups/carnegie_poc/jingze/ppcache.  Entries are keyed by the dump's realpath, size and
mtime, so a rewritten/rewound dump number never returns stale data.  Writes are atomic (tmp + rename)."""
import os, json, hashlib, time
import numpy as np

ROOT = os.environ.get('PPCACHE', '/resnick/groups/carnegie_poc/jingze/ppcache')


# ------------------------------------------------------------------ raw .bin layout
def _parse(path):
    with open(path, 'rb') as fp:
        fp.seek(0, 2); fsize = fp.tell(); fp.seek(0)
        code = fp.readline().split()
        if not code or code[0] != b'Athena':
            raise TypeError('not an Athena .bin: %s' % path)
        ver = code[-1].split(b'=')[-1]
        if ver not in (b'1.1', b'1.2'):
            raise TypeError('unsupported .bin version %s' % ver)
        npc = int(fp.readline().split(b'=')[-1]); ph = {}
        for _ in range(npc - 1):
            k, v = [x.strip() for x in fp.readline().decode().split('=')]
            ph[k] = v
        nvars = int(fp.readline().split(b'=')[-1])
        var_list = [v.decode() for v in fp.readline().split()[1:]]
        hsz = int(fp.readline().split(b'=')[-1])
        header = [ln.decode().split('#')[0].strip() for ln in fp.read(hsz).split(b'\n')]
        header = [ln for ln in header if ln]
        data0 = fp.tell()
        ls, vs = int(ph['size of location']), int(ph['size of variable'])
        coord = int(ph.get('coordinate arrays', '0')) == 1
        ldt = np.float64 if ls == 8 else np.float32; vdt = np.float64 if vs == 8 else np.float32
        # walk the MeshBlock records (headers only)
        recs = []; off = data0
        while off < fsize:
            fp.seek(off)
            idx = np.frombuffer(fp.read(24), np.int32).astype(np.int64)
            lg = np.frombuffer(fp.read(16), np.int32)
            n1, n2, n3 = idx[1] - idx[0] + 1, idx[3] - idx[2] + 1, idx[5] - idx[4] + 1
            hb = 40 + 6 * ls + ((2 * (n1 + n2 + n3) + 3) * ls if coord else 0)
            recs.append((off, idx, lg, int(n1), int(n2), int(n3), hb))
            off += hb + nvars * n1 * n2 * n3 * vs
            if len(recs) == 1 and off < fsize:
                # uniform blocks: the rest follow with the same stride (checked against the file size)
                stride = off - data0
                if (fsize - data0) % stride == 0:
                    nmb = (fsize - data0) // stride
                    fp.seek(data0)
                    raw = np.memmap(path, np.uint8, 'r', offset=data0, shape=(nmb, stride))
                    ok = True
                    for m in range(1, nmb):
                        b = np.asarray(raw[m, :40]).tobytes()
                        i2 = np.frombuffer(b[:24], np.int32).astype(np.int64)
                        if not np.array_equal(i2[1::2] - i2[::2], idx[1::2] - idx[::2]):
                            ok = False; break
                        recs.append((data0 + m * stride, i2, np.frombuffer(b[24:40], np.int32), int(n1), int(n2), int(n3), hb))
                    if ok:
                        break
                    recs = recs[:1]
    return dict(ph=ph, nvars=nvars, var_list=var_list, header=header, recs=recs, ls=ls, vs=vs, coord=coord, ldt=ldt,
                vdt=vdt, fsize=fsize)


def _hval(header, block, key):
    blk = '<none>'
    for ln in header:
        if ln.startswith('<'):
            blk = ln; continue
        body = ln.split('#')[0]
        if '=' in body:
            k, v = body.split('=', 1)
            if blk == block and k.strip() == key:
                return v
    raise KeyError(block + key)


def read_binary(path, vars=None):
    """bin_convert.read_binary with an optional variable subset; identical values and keys."""
    L = _parse(path)
    ph, hdr, recs = L['ph'], L['header'], L['recs']
    want = L['var_list'] if vars is None else [v for v in L['var_list'] if v in vars]
    vi = [L['var_list'].index(v) for v in want]
    mm = np.memmap(path, np.uint8, 'r')
    fd = {'header': hdr, 'time': float(ph['time']), 'cycle': int(ph['cycle']), 'var_names': L['var_list']}
    for k in ('Nx1', 'Nx2', 'Nx3'):
        fd[k] = int(_hval(hdr, '<mesh>', 'n' + k[1:].lower()))
    fd['nvars'] = L['nvars']
    for k in ('x1min', 'x1max', 'x2min', 'x2max', 'x3min', 'x3max'):
        fd[k] = float(_hval(hdr, '<mesh>', k))
    nghost = int(_hval(hdr, '<mesh>', 'nghost'))
    fd['n_mbs'] = len(recs)
    for k in ('nx1', 'nx2', 'nx3'):
        fd[k + '_mb'] = int(_hval(hdr, '<meshblock>', k))
    o0, i0, _, n1, n2, n3, _ = recs[0]
    fd['nx1_out_mb'], fd['nx2_out_mb'], fd['nx3_out_mb'] = n1, n2, n3
    mb_index, mb_logical, mb_geom = [], [], []
    coords = {k: [] for k in ('x1f', 'x1v', 'x2f', 'x2v', 'x3f', 'x3v')}
    mb_data = {v: [] for v in want}
    ls, vs = L['ls'], L['vs']
    for off, idx, lg, n1, n2, n3, hb in recs:
        mb_index.append(idx - nghost); mb_logical.append(np.array(lg))
        mb_geom.append(np.frombuffer(mm[off + 40:off + 40 + 6 * ls].tobytes(), L['ldt']))
        if L['coord']:
            p = off + 40 + 6 * ls
            for d, nout in ((1, n1), (2, n2), (3, n3)):
                for key, n in (('x%df' % d, nout + 1), ('x%dv' % d, nout)):
                    coords[key].append(np.frombuffer(mm[p:p + n * ls].tobytes(), L['ldt'])); p += n * ls
        cell = n1 * n2 * n3
        for v, i in zip(want, vi):
            a = off + hb + i * cell * vs
            mb_data[v].append(np.frombuffer(mm[a:a + cell * vs].tobytes(), L['vdt']).reshape(n3, n2, n1))
    for v in L['var_list']:
        if v not in mb_data:
            mb_data[v] = None if vars is not None else []
    if vars is not None:
        mb_data = {v: mb_data[v] for v in want}
    fd['mb_index'] = np.array(mb_index); fd['mb_logical'] = np.array(mb_logical); fd['mb_geometry'] = np.array(mb_geom)
    fd['mb_data'] = mb_data
    if L['coord']:
        for k, v in coords.items():
            fd['mb_' + k] = np.array(v)
    return fd


# ------------------------------------------------------------------ cache
def _entry(path):
    rp = os.path.realpath(path); st = os.stat(rp)
    run = hashlib.md5(os.path.dirname(rp).encode()).hexdigest()[:12]
    tag = '%s_%d_%d' % (os.path.basename(rp), st.st_size, st.st_mtime_ns)
    d = os.path.join(ROOT, run, tag)
    return d, rp


def _save(fn, a):
    tmp = fn + '.%d.tmp.npy' % os.getpid(); np.save(tmp, a); os.replace(tmp, fn)


def geometry(path):
    d, rp = _entry(path)
    gf = os.path.join(os.path.dirname(d), 'geometry.npz')
    if os.path.exists(gf):
        return dict(np.load(gf, allow_pickle=False))
    h = read_binary(rp, vars=[])
    g = {k: v for k, v in h.items() if isinstance(v, np.ndarray)}
    g.update(time=np.array(h['time']), var_names=np.array(h['var_names']))
    os.makedirs(os.path.dirname(gf), exist_ok=True)
    tmp = gf + '.%d.tmp.npz' % os.getpid(); np.savez(tmp, **g); os.replace(tmp, gf)
    return g


def load(path, fields, cache=True):
    """{field: (nmb, n3, n2, n1) array} (exact dump values) + '_hdr' (dict from read_binary without mb_data)."""
    d, rp = _entry(path)
    out, miss = {}, []
    for f in fields:
        fn = os.path.join(d, f + '.npy')
        if cache and os.path.exists(fn):
            out[f] = np.load(fn, mmap_mode='r')
        else:
            miss.append(f)
    hf = os.path.join(d, 'hdr.json')
    if miss or not os.path.exists(hf):
        h = read_binary(rp, vars=miss)
        for f in miss:
            out[f] = np.stack(h['mb_data'][f])
        hdr = {k: v for k, v in h.items() if k != 'mb_data'}
        if cache:
            os.makedirs(d, exist_ok=True)
            for f in miss:
                _save(os.path.join(d, f + '.npy'), out[f])
            js = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in hdr.items()}
            tmp = hf + '.%d.tmp' % os.getpid(); json.dump(js, open(tmp, 'w')); os.replace(tmp, hf)
    else:
        hdr = json.load(open(hf))
    if cache:
        try:
            os.utime(d)
        except OSError:
            pass
    out['_hdr'] = hdr
    return out


def derived(path, name, version, fn, deps, dtype=None):
    d, rp = _entry(path)
    key = '%s__%s' % (name, hashlib.md5(str(version).encode()).hexdigest()[:8])
    mf = os.path.join(d, key + '.json')
    if os.path.exists(mf):
        names = json.load(open(mf))
        r = {k: np.load(os.path.join(d, '%s.%s.npy' % (key, k)), mmap_mode='r') for k in names}
        return r if len(names) != 1 or names[0] != '_' else r['_']
    src = load(rp, deps)
    r = fn({k: v for k, v in src.items() if k != '_hdr'})
    rd = r if isinstance(r, dict) else {'_': r}
    os.makedirs(d, exist_ok=True)
    for k, a in rd.items():
        a = np.asarray(a)
        _save(os.path.join(d, '%s.%s.npy' % (key, k)), a.astype(dtype) if dtype else a)
    tmp = mf + '.%d.tmp' % os.getpid(); json.dump(list(rd), open(tmp, 'w')); os.replace(tmp, mf)
    return r


def prune(max_gb=300.0):
    ents = []
    for run in os.listdir(ROOT) if os.path.isdir(ROOT) else []:
        rd = os.path.join(ROOT, run)
        for t in os.listdir(rd):
            p = os.path.join(rd, t)
            if os.path.isdir(p):
                sz = sum(os.path.getsize(os.path.join(p, f)) for f in os.listdir(p))
                ents.append((os.path.getmtime(p), sz, p))
    ents.sort(reverse=True); tot = 0; gone = 0
    import shutil
    for mt, sz, p in ents:
        tot += sz
        if tot > max_gb * 1e9:
            shutil.rmtree(p, ignore_errors=True); gone += 1
    return gone


# ------------------------------------------------------------------ ready-made derived fields
def dhj_Tp(path, eosf='/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt'):
    """(T, p) of a dhj dump from the run's EOS table, exactly as dhjcs.EOS(eosf).invert(dens, eint) (float64, cached)."""
    def fn(F):
        import sys
        sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rochedhj_ana'); import dhjcs
        T, p = dhjcs.EOS(eosf).invert(np.asarray(F['dens']), np.asarray(F['eint']))
        return {'T': T, 'p': p}
    st = os.stat(eosf)
    return derived(path, 'dhj_Tp', '%s:%d:%d:v1' % (os.path.realpath(eosf), st.st_size, st.st_mtime_ns), fn, ['dens', 'eint'])
