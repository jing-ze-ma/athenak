"""Fast partial reader for cyclediag: reads only selected variables' (k,j,:) rows."""
import numpy as np


def header(path):
    with open(path, 'rb') as f:
        h = f.readline()
        n = f.readline()
        base = f.tell()
    hdr = h.decode().strip()
    names = n.decode().strip().split()
    out = {}
    for tok in hdr.split()[2:]:
        key, _, val = tok.partition('=')
        try:
            out[key] = int(val)
        except ValueError:
            out[key] = float(val)
    out['names'] = names
    out['base'] = base
    return out


FACE = ('rad_w', 'rad_tauf')
COL = ('rt_icut',)


def offsets(h):
    n1, n2, n3, nf = h['n1'], h['n2'], h['n3'], h['nface']
    off = {}
    o = 0
    for nm in h['names']:
        if nm in COL:
            sz = n3*n2
        elif nm in FACE:
            sz = n3*n2*nf
        else:
            sz = n3*n2*n1
        off[nm] = (o, sz)
        o += sz
    return off


def read_rows(path, vars, k, j):
    h = header(path)
    off = offsets(h)
    n1, n2 = h['n1'], h['n2']
    res = {'cycle': h['cycle'], 'time': h['time'], 'dt': h['dt'], 'is': h['is'],
           'ie': h['ie'],
           'js': h['js'], 'je': h['je'], 'ks': h['ks'], 'ke': h['ke'], 'n1': n1, 'n2': n2,
           'n3': h['n3']}
    with open(path, 'rb') as f:
        for v in vars:
            o, sz = off[v]
            nn = h['nface'] if v in FACE else n1
            start = o+(k*n2+j)*nn
            f.seek(h['base']+start*8)
            res[v] = np.frombuffer(f.read(nn*8), dtype='<f8')
    return res


def read_full(path, vars):
    h = header(path)
    off = offsets(h)
    n1, n2, n3 = h['n1'], h['n2'], h['n3']
    res = {'cycle': h['cycle'], 'time': h['time'], 'dt': h['dt'], 'is': h['is'],
           'ie': h['ie'],
           'js': h['js'], 'je': h['je'], 'ks': h['ks'], 'ke': h['ke']}
    with open(path, 'rb') as f:
        for v in vars:
            o, sz = off[v]
            f.seek(h['base']+o*8)
            res[v] = np.frombuffer(f.read(sz*8), dtype='<f8').reshape((n3, n2, -1))
    return res
