"""Cubed-sphere grid imprint of rho/<rho> - 1 in star3d_cs dumps.

usage: python3 cshero_imprint.py DUMP.hydro_w.NNNNN.bin [...]   (numpy only)
Readers copied from the star3d notebook (athenak-vet-gallery star3d-cs).
Prints per radius bin: equatorial m = 4 / m = 8 amplitude (panels 0,1,2,4 at eta = 0) and
the panel-symmetric rms (drho averaged over the 6 panels at each lateral index, minus its
solid-angle mean): the part of drho that has the cube's symmetry.
"""
import sys
import numpy as np
import os
import re


def read_bin(path):
    """One .bin dump: dict with 'time', 'cycle', 'vars', 'deck' (text) and 'blocks', a list
    (in MeshBlock id order) of dicts: lx (logical location), x1v, x2v, x3v (cell centres),
    and one (nx3, nx2, nx1) array per variable."""
    with open(path, 'rb') as fp:
        size = os.fstat(fp.fileno()).st_size
        assert fp.readline().split()[0] == b'Athena', path
        npre = int(fp.readline().split(b'=')[-1])
        par = {}
        for _ in range(npre - 1):
            k, v = [x.strip() for x in fp.readline().decode().split('=')]
            par[k] = v
        coords = int(par.get('coordinate arrays', '0')) == 1
        lsz, vsz = int(par['size of location']), int(par['size of variable'])
        fp.readline()
        names = fp.readline().decode().split()[1:]
        deck = fp.read(int(fp.readline().split(b'=')[-1])).decode(errors='replace')
        ldt = np.float64 if lsz == 8 else np.float32
        vdt = np.float64 if vsz == 8 else np.float32
        blocks = []
        while fp.tell() < size:
            ind = np.frombuffer(fp.read(24), dtype=np.int32)
            lx = np.frombuffer(fp.read(16), dtype=np.int32)
            ext = np.frombuffer(fp.read(6*lsz), dtype=ldt)
            n1, n2, n3 = ind[1]-ind[0]+1, ind[3]-ind[2]+1, ind[5]-ind[4]+1
            b = {'lx': lx[:3].copy(), 'ext': ext.copy()}
            if coords:
                for d, n in (('1', n1), ('2', n2), ('3', n3)):
                    b['x%sf' % d] = np.fromfile(fp, dtype=ldt, count=n + 1)
                    b['x%sv' % d] = np.fromfile(fp, dtype=ldt, count=n)
            else:
                for d, n, lo, hi in (('1', n1, ext[0], ext[1]), ('2', n2, ext[2], ext[3]),
                                     ('3', n3, ext[4], ext[5])):
                    b['x%sf' % d] = lo + np.arange(n + 1)*(hi - lo)/n
                    b['x%sv' % d] = 0.5*(b['x%sf' % d][1:] + b['x%sf' % d][:-1])
            a = np.fromfile(fp, dtype=vdt, count=n1*n2*n3*len(names))
            a = a.reshape(len(names), n3, n2, n1)
            for i, n in enumerate(names):
                b[n] = a[i]
            blocks.append(b)
    return {'time': float(par['time']), 'cycle': int(par['cycle']), 'vars': names,
            'deck': deck, 'blocks': blocks}


def assemble(dump, var):
    """Global array of var: (nx3, nx2, nx1) on the sp wedge, (6, nx3, nx2, nx1) on the cubed
    sphere (MeshBlocks are written panel by panel, so panel = block index // blocks per
    panel), plus the radial cell faces x1f (one MeshBlock spans the whole radius)."""
    bl = dump['blocks']
    n3, n2, n1 = bl[0][var].shape
    cs = re.search(r'^\s*use_cubed_sphere\s*=\s*(true|1)', dump['deck'], re.M) is not None
    npan = 6 if cs else 1
    per = len(bl)//npan
    L2 = max(b['lx'][1] for b in bl) + 1
    L3 = max(b['lx'][2] for b in bl) + 1
    out = np.zeros((npan, L3*n3, L2*n2, n1))
    for m, b in enumerate(bl):
        p = m//per
        j0, k0 = b['lx'][1]*n2, b['lx'][2]*n3
        out[p, k0:k0+n3, j0:j0+n2, :] = b[var]
    return (out if cs else out[0]), bl[0]['x1f'], cs

RSUN = 6.957e10


def solid_angle(n):
    """Solid angle of the n x n equiangular cells of one panel, (eta, xi)."""
    e = np.tan(np.linspace(-np.pi/4, np.pi/4, n + 1))
    X, Y = np.meshgrid(e, e)
    f = np.arctan2(X*Y, np.sqrt(1 + X*X + Y*Y))   # corner function of the gnomonic map
    return f[1:, 1:] - f[1:, :-1] - f[:-1, 1:] + f[:-1, :-1]


def analyse(path):
    d = read_bin(path)
    rho, rf, cs = assemble(d, 'dens')
    n = rho.shape[1]
    w = solid_angle(n)[None, :, :, None]
    mean = (rho*w).sum(axis=(0, 1, 2))/(6*w).sum(axis=(0, 1, 2))
    drho = rho/mean - 1
    sym = drho.mean(axis=0)                                    # (eta, xi, r)
    sym = sym - (sym*w[0]).sum(axis=(0, 1))/w[0].sum()
    rms_sym = np.sqrt((sym**2*w[0]).sum(axis=(0, 1))/w[0].sum())
    rms_all = np.sqrt((drho**2*w).sum(axis=(0, 1, 2))/(6*w).sum())
    ring = np.concatenate([drho[p, n//2-1:n//2+1].mean(0) for p in (0, 1, 2, 4)], axis=0)
    c = np.fft.rfft(ring, axis=0)/ring.shape[0]
    a4, a8 = 2*np.abs(c[4]), 2*np.abs(c[8])
    r = 0.5*(rf[1:] + rf[:-1])/RSUN
    return d['cycle'], d['time'], r, a4, a8, rms_sym, rms_all


if __name__ == '__main__':
    edges = [0, 100, 200, 400, 600, 800, 850, 900, 1000, 2000]
    for p in sys.argv[1:]:
        cyc, t, r, a4, a8, rs, ra = analyse(p)
        print('%s cycle %d t = %.3f d' % (p, cyc, t/86400))
        print('  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all')
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = (r >= lo) & (r < hi)
            if s.any():
                print('  %4d-%-5d  %.3e  %.3e  %.3e  %.3e' % (lo, hi, a4[s].max(), a8[s].max(),
                      rs[s].max(), ra[s].max()))
        np.savez(p + '.imprint.npz', r=r, a4=a4, a8=a8, rms_sym=rs, rms_all=ra, cycle=cyc, time=t)
