"""libin.py: read AthenaK .bin outputs (version 1.1/1.2, float or double) and assemble the
MeshBlocks into GLOBAL arrays indexed by (panel, k, j, i), so that two runs with different
MeshBlock layouts / rank counts can be compared cell by cell.

Only uniform (single-level) meshes.  Cubed sphere: MeshBlock gids are panel-major (one
tree per panel, mesh/build_tree.cpp), so panel = gid // (nmb/6) and lx2/lx3 are
panel-local.  The file stores MeshBlocks in gid order (rank offsets by gid)."""
import os
import numpy as np


def read_bin(fn):
    with open(fn, 'rb') as f:
        raw = f.read()
    pos = 0
    def line():
        nonlocal pos
        e = raw.index(b'\n', pos)
        s = raw[pos:e].decode()
        pos = e + 1
        return s
    ver = line()
    hdr = {}
    while True:
        s = line().strip()
        if s.startswith('variables:'):
            hdr['variables'] = s.split(':', 1)[1].split()
            break
        k, v = s.split('=')
        hdr[k.strip()] = v.strip()
    s = line().strip()
    hoff = int(s.split('=')[1])
    pin_txt = raw[pos:pos + hoff].decode(errors='replace')
    pos += hoff
    lsz = int(hdr['size of location'])
    vsz = int(hdr['size of variable'])
    coord = hdr.get('coordinate arrays', '0') == '1'
    nvar = int(hdr['number of variables'])
    ldt = np.float64 if lsz == 8 else np.float32
    vdt = np.float64 if vsz == 8 else np.float32
    blocks = []
    while pos < len(raw):
        ii = np.frombuffer(raw, np.int32, 10, pos)
        pos += 40
        ois, oie, ojs, oje, oks, oke, lx1, lx2, lx3, lev = [int(x) for x in ii]
        n1, n2, n3 = oie - ois + 1, oje - ojs + 1, oke - oks + 1
        pos += 6 * lsz
        if coord:
            pos += ((2 * n1 + 1) + (2 * n2 + 1) + (2 * n3 + 1)) * lsz
        cells = n1 * n2 * n3
        a = np.frombuffer(raw, vdt, nvar * cells, pos).reshape(nvar, n3, n2, n1)
        pos += nvar * cells * vsz
        blocks.append(dict(l=(lx1, lx2, lx3), lev=lev, n=(n1, n2, n3), data=a))
    pin = parse_pin(pin_txt)
    return dict(time=float(hdr['time']), cycle=int(hdr['cycle']), vars=hdr['variables'],
                blocks=blocks, pin=pin, vsz=vsz)


def parse_pin(txt):
    d, blk = {}, None
    for s in txt.splitlines():
        s = s.split('#', 1)[0].strip()
        if not s:
            continue
        if s.startswith('<'):
            blk = s.strip('<>').strip()
            d.setdefault(blk, {})
        elif '=' in s and blk is not None:
            k, v = s.split('=', 1)
            d[blk][k.strip()] = v.strip()
    return d


def assemble(b):
    """-> (G[nvar, npanel, N3, N2, N1], layout) ; layout = dict(nmb, n=(n1,n2,n3), npanel)"""
    blks = b['blocks']
    if any(x['lev'] != 0 for x in blks):
        raise RuntimeError('refined mesh: libin.assemble handles single-level meshes only')
    m = b['pin'].get('mesh', {})
    cs = m.get('use_cubed_sphere', 'false').lower() in ('true', '1')
    npan = 6 if cs else 1
    nmb = len(blks)
    if nmb % npan:
        raise RuntimeError('nmb %d not a multiple of %d panels' % (nmb, npan))
    per = nmb // npan
    n1, n2, n3 = blks[0]['n']
    N1 = int(m['nx1']); N2 = int(m['nx2']); N3 = int(m['nx3'])
    # output may be a slice / subset: size the global array from the blocks
    mx = [0, 0, 0]
    for x in blks:
        for d in range(3):
            mx[d] = max(mx[d], (x['l'][d] + 1) * x['n'][d])
    G = np.full((blks[0]['data'].shape[0], npan, mx[2], mx[1], mx[0]), np.nan)
    for g, x in enumerate(blks):
        p = g // per
        l1, l2, l3 = x['l']
        G[:, p, l3 * n3:(l3 + 1) * n3, l2 * n2:(l2 + 1) * n2, l1 * n1:(l1 + 1) * n1] = x['data']
    if np.isnan(G).any():
        raise RuntimeError('assembly left holes: layout not recognised')
    return G, dict(nmb=nmb, n=(n1, n2, n3), npanel=npan, N=(N1, N2, N3))


def files_by_cycle(d):
    """map cycle -> list of .bin files (one per output variable) in run dir d/bin"""
    out = {}
    bd = os.path.join(d, 'bin')
    if not os.path.isdir(bd):
        return out
    for fn in sorted(os.listdir(bd)):
        if not fn.endswith('.bin'):
            continue
        p = os.path.join(bd, fn)
        with open(p, 'rb') as f:
            head = f.read(400).decode(errors='replace')
        cyc = int(head.split('cycle=')[1].split()[0])
        out.setdefault(cyc, []).append(p)
    return out


def load_cycle(d, cyc):
    """all variables of run d at cycle cyc -> (dict name -> global array, layout)"""
    res, lay = {}, None
    for p in files_by_cycle(d).get(cyc, []):
        b = read_bin(p)
        G, lay = assemble(b)
        for v, nm in enumerate(b['vars']):
            res[nm] = G[v]
    return res, lay
