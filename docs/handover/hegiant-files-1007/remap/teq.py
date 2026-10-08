#!/usr/bin/env python3
"""Gas-radiation temperature mismatch |T_rad/T_gas - 1| by radial zone on 64 sample columns: T_gas from the code's
EOS table dump (e_int = E - KE - rho Phi with the he_gm_column Phi), T_rad = (E_rad/a)^(1/4).
usage: teq.py A.rst [B.rst ...]"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_giant_1006/remap')
import he_remap_rst_w256 as H
from phi_code import phi_cells
RS, AR = 6.957e10, 7.5657332503e-15
fn = '/viper/ptmp2/jinma/he_giant_1006/eos/eos_dump_y0.9265_a1.77.txt'
with open(fn) as fh:
    ln = [fh.readline() for _ in range(3)]
nx, ny, x0, dx, y0, dy = [float(v) for v in ln[2].split()[1:]]
nx, ny = int(nx), int(ny)
d = np.loadtxt(fn, comments='#')
LE = d[:, 0].reshape(ny, nx)
lT = y0 + dy * np.arange(ny)


def tgas(rho, eint):
    lr = np.log10(rho)
    i = np.clip(((lr - x0) / dx).astype(int), 0, nx - 2)
    w = (lr - x0) / dx - i
    out = np.empty(rho.size)
    target = np.log10(eint / rho)
    for k in range(rho.size):
        col = (1 - w[k]) * LE[:, i[k]] + w[k] * LE[:, i[k] + 1]
        out[k] = 10**np.interp(target[k], col, lT)
    return out


for f in sys.argv[1:]:
    r = H.read_rst(f)
    G = H.gather(r)
    ng, n1 = r.ng, r.bind[1]
    p = [float(H.get_param(r.text, 'mesh', k) or 0.0) for k in H.SKEYS]
    fa = H.faces_of(p, n1, ng, r.msize[0], r.msize[3])
    x1v, ph = phi_cells(fa, ng, H.get_param(r.text, 'problem', 'he_ic_file'),
                        float(H.get_param(r.text, 'problem', 'he_gm')), int(H.get_param(r.text, 'problem', 'he_nfine')))
    a = slice(ng, ng + n1)
    u = G['hyd'][:, ::16, ::16, a]
    E = G['m1'][0][::16, ::16, a]
    ke = 0.5 * (u[1]**2 + u[2]**2 + u[3]**2) / u[0]
    ei = u[4] - ke - u[0] * ph[a]
    Tg = tgas(u[0].ravel(), np.maximum(ei.ravel(), 1e-300)).reshape(u[0].shape)
    Tr = (E / AR)**0.25
    dd = np.abs(Tr / Tg - 1)
    x = x1v[a] / RS
    s = f.split('/')[-2] + ':'
    for lo, hi in [(3, 30), (30, 55), (55, 62), (62, 65), (65, 70), (70, 200.1)]:
        m = (x >= lo) & (x < hi)
        s += ' | %g-%g med %.1e max %.1e' % (lo, hi, np.median(dd[..., m]), dd[..., m].max())
    print(s)
