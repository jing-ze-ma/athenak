"""Gas temperature and gas pressure from (rho, eint) with the run's OWN general-EOS table (Saha H/He, X 0.70 Y 0.28,
ionization on, no radiation in the EOS): the <hydro>/eos_table_dump of the production EOS keys (EOSTBL, md5 EOSMD5).
Table: log10(e/rho), log10(p/rho) [cgs] on (log rho 901 x log T 1026, T slowest).
T: per log-rho column, T(log e) on a uniform log-e grid (monotone interp), then bilinear in (log rho, log e).
P_gas: bilinear in (log rho, log T) of the table's log p/rho at that T. Inputs must be float64 (float32 input rounds the
clip ne-1.000001 to ne-1 -> IndexError). Parsed table cached in $OUT/data/eos_cache_<md5>.npz."""
import os
import numpy as np
from common import EOSTBL, EOSMD5, OUT, RHO_U, E_U, md5


class TGas:
    def __init__(self, fn=EOSTBL, ne=4000):
        if not os.path.exists(fn):
            raise SystemExit(f'EOS table {fn} not found: set EOSTBL (see mk_eostbl.sh to regenerate it, md5 {EOSMD5})')
        m = md5(fn)
        if m != EOSMD5:
            print(f'WARNING EOS table md5 {m} != expected {EOSMD5} ({fn})', flush=True)
        self.md5 = m
        cache = os.path.join(OUT, 'data', f'eos_cache_{m[:8]}_{ne}.npz')
        if os.path.exists(cache):
            c = np.load(cache)
            for k in c.files:
                setattr(self, k, c[k] if c[k].ndim else c[k].item())
            return
        hdr = [ln for ln in open(fn) if ln.startswith('#')][:6]
        nx, ny, x0, dx, y0, dy = [float(v) for v in hdr[2][1:].split()]
        nx, ny = int(nx), int(ny)
        d = np.loadtxt(fn).reshape(ny, nx, 2)          # [T, rho, (log e/rho, log p/rho)]
        self.x0, self.dx, self.nx, self.y0, self.dy, self.ny = x0, dx, nx, y0, dy, ny
        ylog = y0 + dy * np.arange(ny)
        le = np.maximum.accumulate(d[..., 0], axis=0)  # guard monotonicity
        self.lp = d[..., 1].copy()
        self.e0, self.e1 = float(le.min()), float(le.max())
        self.ne = ne
        self.eg = np.linspace(self.e0, self.e1, ne)
        self.de = self.eg[1] - self.eg[0]
        inv = np.empty((nx, ne), np.float64)
        for i in range(nx):
            inv[i] = np.interp(self.eg, le[:, i], ylog)
        self.inv = inv
        self.ylo, self.yhi = ylog[0], ylog[-1]
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        np.savez(cache + '.tmp.npz', **{k: getattr(self, k) for k in
                                        ('x0', 'dx', 'nx', 'y0', 'dy', 'ny', 'e0', 'e1', 'ne', 'eg', 'de', 'inv', 'ylo', 'yhi', 'lp')})
        os.replace(cache + '.tmp.npz', cache)

    def logT(self, rho_code, eint_code):
        x = np.log10(np.maximum(rho_code, 1e-40) * RHO_U)
        e = np.log10(np.maximum(eint_code / np.maximum(rho_code, 1e-40), 1e-40) * E_U)
        fx = np.clip((x - self.x0) / self.dx, 0, self.nx - 1.000001)
        fe = np.clip((e - self.e0) / self.de, 0, self.ne - 1.000001)
        i = fx.astype(np.int64); j = fe.astype(np.int64)
        wx = fx - i; we = fe - j
        Iv = self.inv
        return ((1 - wx) * ((1 - we) * Iv[i, j] + we * Iv[i, j + 1]) + wx * ((1 - we) * Iv[i + 1, j] + we * Iv[i + 1, j + 1]))

    def logPgas(self, rho_code, logT):
        """log10 P_gas in code units (rho_code (km/s)^2) = log10 rho_code + log10(p/rho)_cgs - 10."""
        x = np.log10(np.maximum(rho_code, 1e-40) * RHO_U)
        fx = np.clip((x - self.x0) / self.dx, 0, self.nx - 1.000001)
        fy = np.clip((logT - self.y0) / self.dy, 0, self.ny - 1.000001)
        i = fx.astype(np.int64); j = fy.astype(np.int64)
        wx = fx - i; wy = fy - j
        L = self.lp
        lpr = (1 - wy) * ((1 - wx) * L[j, i] + wx * L[j, i + 1]) + wy * ((1 - wx) * L[j + 1, i] + wx * L[j + 1, i + 1])
        return np.log10(np.maximum(rho_code, 1e-40)) + lpr - np.log10(E_U)
