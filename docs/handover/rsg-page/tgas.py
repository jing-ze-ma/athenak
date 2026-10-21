#!/usr/bin/env python3
"""T_gas(rho, eint) for the RSG E deck, from the code's OWN general-EOS table dump (<hydro>/eos_table_dump).

eos_dump_rsgE.txt = dump of eos = general, general_eos = table, eos_ionization true, eos_h2 false, eos_metal_ionization false,
eos_radiation false, X 0.653848, Y 0.326395, a_metal 16, logd -21..-1, logT 3.3..6.8 (Raven ic_star_cache eos_85908b8dfca6).
Units: the RSG deck is cgs (<units> all 1), so hydro_w dens [g/cm^3] and eint [erg/cm^3 per volume] go in directly.
Method: per log-rho column, invert the monotone log10(e/rho)(log T) for log T, then interpolate linearly in log rho
(the code uses bicubic Hermite on the same nodes; the difference is far below plotting accuracy). T is clipped to the
table range; apply tfloor 2000 K yourself if you want the code's floor.
Usage: from tgas import TGas; tg = TGas('eos_dump_rsgE.txt'); T = tg(rho, eint)
"""
import numpy as np


class TGas:
    def __init__(self, fn):
        with open(fn) as f:
            f.readline(); f.readline()
            nx, ny, x0, dx, y0, dy = f.readline()[1:].split()
        self.nx, self.ny = int(nx), int(ny)
        self.x0, self.dx, self.y0, self.dy = float(x0), float(dx), float(y0), float(dy)
        d = np.loadtxt(fn)
        self.le = d[:, 0].reshape(self.ny, self.nx)          # log10(e/rho) [T row, rho col]
        self.lt = self.y0 + self.dy * np.arange(self.ny)

    def _col(self, j, les):
        c = self.le[:, j]
        return np.interp(les, c, self.lt)                     # monotone in T

    def __call__(self, rho, eint):
        rho = np.asarray(rho, dtype=float); eint = np.asarray(eint, dtype=float)
        lr = np.clip(np.log10(rho), self.x0, self.x0 + self.dx * (self.nx - 1))
        les = np.log10(np.maximum(eint / rho, 1e-300))
        fx = (lr - self.x0) / self.dx
        j = np.clip(np.floor(fx).astype(int), 0, self.nx - 2); w = fx - j
        out = np.empty_like(lr)
        flat_j, flat_w, flat_e = j.ravel(), w.ravel(), les.ravel()
        res = out.ravel()
        for jj in np.unique(flat_j):
            m = flat_j == jj
            t0 = self._col(jj, flat_e[m]); t1 = self._col(jj + 1, flat_e[m])
            res[m] = (1 - flat_w[m]) * t0 + flat_w[m] * t1
        return 10.0 ** res.reshape(lr.shape)


if __name__ == '__main__':
    import sys
    tg = TGas(sys.argv[1] if len(sys.argv) > 1 else 'eos_dump_rsgE.txt')
    # self-check on table nodes: (rho, T) -> e from the table -> back to T
    for lr, lt in [(-9.0, 4.0), (-12.0, 3.7), (-6.0, 5.5), (-15.0, 3.5)]:
        j = int(round((lr - tg.x0) / tg.dx)); i = int(round((lt - tg.y0) / tg.dy))
        e = 10 ** tg.le[i, j] * 10 ** lr
        print('rho 1e%.1f T %.0f K -> eint %.6e erg/cm3 -> T %.1f K' % (lr, 10 ** lt, e, tg(10 ** lr, e)))
