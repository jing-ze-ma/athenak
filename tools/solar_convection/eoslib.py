"""Offline evaluation of AthenaK's general EOS, for analysing solar_convection dumps.

The solver has no temperature output for Newtonian hydro, and the tabulated EOS lives in a
device array, so the only faithful way to turn the (rho, e) pairs in a .bin dump into
T / mu / Gamma_1 / grad_ad is to re-evaluate the SAME composition model offline.
tools/eos_dump.cpp includes the code's own eos_composition.hpp and writes a
(log10 rho, log10 T) grid; everything here interpolates that grid.

Both EOS choices are exposed through one interface so the analysis script does not branch:

    eos = IdealEOS(Rgas=1.38e8, gamma=5/3)      # the control run
    eos = TableEOS('eos_grid.bin')              # general_eos = table

Every method takes cgs density and cgs SPECIFIC internal energy e/rho and returns cgs.
"""

import os

import numpy as np

# ---------------------------------------------------------------------------------------
# ideal gas: composition is fixed by Rgas, so everything is a closed form


class IdealEOS:
    name = 'ideal'

    def __init__(self, Rgas=1.38e8, gamma=5.0/3.0):
        self.Rgas = Rgas
        self.gamma = gamma
        # k_B/(m_u Rgas) is the mean molecular weight the pgen's Rgas implies
        self.mu_fixed = 1.380649e-16/(1.66053906660e-24*Rgas)

    def temperature(self, rho, espec):
        return (self.gamma - 1.0)*espec/self.Rgas

    def pressure(self, rho, espec):
        return (self.gamma - 1.0)*rho*espec

    def mu(self, rho, espec):
        return np.full_like(np.asarray(rho, dtype=float), self.mu_fixed)

    def gamma1(self, rho, espec):
        return np.full_like(np.asarray(rho, dtype=float), self.gamma)

    def grad_ad(self, rho, espec):
        return np.full_like(np.asarray(rho, dtype=float),
                            (self.gamma - 1.0)/self.gamma)


# ---------------------------------------------------------------------------------------
# tabulated general EOS


def _bilinear(field, ix, iy, wx, wy):
    return ((1.0-wx)*(1.0-wy)*field[ix, iy] + wx*(1.0-wy)*field[ix+1, iy] +
            (1.0-wx)*wy*field[ix, iy+1] + wx*wy*field[ix+1, iy+1])


class TableEOS:
    name = 'table'

    def __init__(self, path):
        raw = np.fromfile(path, dtype=np.float64)
        nd, nt = int(raw[0]), int(raw[1])
        self.logd_min, self.dlogd, self.logt_min, self.dlogt = raw[2:6]
        self.nd, self.nt = nd, nt
        body = raw[6:].reshape(nd, nt, 5)
        self.log_e = np.ascontiguousarray(body[:, :, 0])   # log10 e_spec [erg/g]
        self.log_p = np.ascontiguousarray(body[:, :, 1])   # log10 (p/rho) [erg/g]
        self._mu = np.ascontiguousarray(body[:, :, 2])
        self.xh2 = np.ascontiguousarray(body[:, :, 3])
        self.xhii = np.ascontiguousarray(body[:, :, 4])
        self.logd = self.logd_min + self.dlogd*np.arange(nd)
        self.logt = self.logt_min + self.dlogt*np.arange(nt)

        # Thermodynamic derivatives, from the log-log grid.  ln10 factors cancel because
        # every derivative below is logarithmic in both variables.
        ln10 = np.log(10.0)
        # chi_rho = (dln p/dln rho)_T,  chi_T = (dln p/dln T)_rho.  p = rho * p_spec, so
        # dln p/dln rho = 1 + dln p_spec/dln rho.
        dlp_dld = np.gradient(self.log_p, self.dlogd, axis=0)
        dlp_dlt = np.gradient(self.log_p, self.dlogt, axis=1)
        dle_dlt = np.gradient(self.log_e, self.dlogt, axis=1)
        self.chi_rho = 1.0 + dlp_dld
        self.chi_T = dlp_dlt
        # cv = (de/dT)_rho = (e/T) dln e/dln T
        e = 10.0**self.log_e
        t2d = 10.0**self.logt[None, :]
        self.cv = e/t2d*dle_dlt
        # p/(rho T cv) = p_spec/(T cv)
        p_spec = 10.0**self.log_p
        aux = p_spec/(t2d*self.cv)
        self.g1 = self.chi_rho + aux*self.chi_T**2
        self.g3m1 = aux*self.chi_T
        self.gradad = self.g3m1/self.g1
        del ln10

    # -- grid lookup ---------------------------------------------------------------------
    def _weights(self, logd, logt):
        fx = np.clip((logd - self.logd_min)/self.dlogd, 0.0, self.nd - 1.0000001)
        fy = np.clip((logt - self.logt_min)/self.dlogt, 0.0, self.nt - 1.0000001)
        ix = fx.astype(np.int64)
        iy = fy.astype(np.int64)
        return ix, iy, fx - ix, fy - iy

    def _eval(self, field, logd, logt):
        ix, iy, wx, wy = self._weights(logd, logt)
        return _bilinear(field, ix, iy, wx, wy)

    # -- inversion -----------------------------------------------------------------------
    def temperature(self, rho, espec):
        """T [K] from cgs density and specific internal energy, by bisection on log10 T.

        e is monotone increasing in T at fixed rho for this model (translational +
        dissociation + ionization all increase with T), so plain bisection is safe.
        """
        logd = np.log10(np.asarray(rho, dtype=float))
        target = np.log10(np.asarray(espec, dtype=float))
        lo = np.full(target.shape, self.logt_min)
        hi = np.full(target.shape, self.logt[-1])
        for _ in range(60):
            mid = 0.5*(lo + hi)
            f = self._eval(self.log_e, logd, mid)
            too_cold = f < target
            lo = np.where(too_cold, mid, lo)
            hi = np.where(too_cold, hi, mid)
        return 10.0**(0.5*(lo + hi))

    # -- accessors, all taking (rho, espec) so the interface matches IdealEOS ------------
    def _ld_lt(self, rho, espec):
        return (np.log10(np.asarray(rho, dtype=float)),
                np.log10(self.temperature(rho, espec)))

    def pressure(self, rho, espec):
        ld, lt = self._ld_lt(rho, espec)
        return np.asarray(rho, dtype=float)*10.0**self._eval(self.log_p, ld, lt)

    def mu(self, rho, espec):
        ld, lt = self._ld_lt(rho, espec)
        return self._eval(self._mu, ld, lt)

    def gamma1(self, rho, espec):
        ld, lt = self._ld_lt(rho, espec)
        return self._eval(self.g1, ld, lt)

    def grad_ad(self, rho, espec):
        ld, lt = self._ld_lt(rho, espec)
        return self._eval(self.gradad, ld, lt)

    def ion_frac(self, rho, espec):
        """(H2 fraction, H+ fraction) of the hydrogen nuclei."""
        ld, lt = self._ld_lt(rho, espec)
        return self._eval(self.xh2, ld, lt), self._eval(self.xhii, ld, lt)


# ---------------------------------------------------------------------------------------
# the pgen's own opacity, so tau is computed the same way the run computes it


def get_kapr(rho, T):
    """Composite H- / Kramers / electron-scattering opacity with the 1e-2 floor.

    Mirrors get_kapr() in src/pgen/solar_convection.cpp exactly (C = 20).
    """
    rho = np.asarray(rho, dtype=float)
    T = np.asarray(T, dtype=float)
    T4 = T/1.0e4
    kap = 2.0e1*np.sqrt(rho/1.0e-6)*T4**9          # H-
    kapd = 3.0e3*(rho/1.0e-5)/np.sqrt((T/1.0e5)**7)  # Kramers
    kap = np.minimum(kap, kapd)
    kap = kap + np.where(T < 1.0e4, 0.0, 0.34)     # electron scattering
    return np.maximum(kap, 1.0e-2)


def eos_for_run(run_dir, grid, athinput='sun.athinput'):
    """The EOS a run actually used, read from its own input file.

    Guessing from the run's NAME is how a tabulated-EOS run gets silently analysed with
    ideal-gas formulas -- 'sponge08' contains neither 'table' nor 'ideal'. Parse
    <hydro>/general_eos instead, and fail loudly rather than fall back to a guess.
    """
    path = os.path.join(run_dir, athinput)
    if not os.path.exists(path):
        raise RuntimeError('%s not found: cannot tell which EOS %s used'
                           % (path, run_dir))
    block, general, kind = None, None, None
    with open(path) as fh:
        for line in fh:
            line = line.split('#')[0].strip()
            if line.startswith('<') and line.endswith('>'):
                block = line[1:-1].strip()
            elif '=' in line and block == 'hydro':
                key, val = (s.strip() for s in line.split('=', 1))
                if key == 'eos':
                    general = val
                elif key == 'general_eos':
                    kind = val
    if general == 'general':
        if kind != 'table':
            raise RuntimeError('%s uses general_eos = %r, which eoslib cannot evaluate'
                               % (run_dir, kind))
        return TableEOS(grid)
    return IdealEOS()
