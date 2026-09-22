"""Q1/Q3: reproduce the per-cell MHD CFL limit of mhd_newdt.cpp from the last dump,
and locate the floored cells.  Read-only."""
import sys
import glob
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert                       # noqa: E402
import dhjcs                             # noqa: E402
import common as C                       # noqa: E402

EOS = C.EOS2()


def load(path):
    raw = bin_convert.read_binary(path)
    mb = raw['mb_data']
    d = {k: np.asarray(mb[k], dtype=np.float64) for k in raw['var_names']}
    g = C.cs_metric(raw)
    geo = dhjcs.cs_geometry(raw)
    T, p, g1 = EOS.invert(d['dens'], d['eint'])
    return raw, d, g, geo, T, p, g1


def cell_dt(d, g, T, p, g1):
    """dt per cell exactly as mhd_newdt.cpp does it (before cfl_no)."""
    rho = d['dens']
    b1, b2, b3 = d['bcc1'], d['bcc2'], d['bcc3']
    v1, v2, v3 = d['velx'], d['vely'], d['velz']
    sn = g['sin_cell'][..., None]
    cf1 = C.fast_speed(rho, p, g1, b1, b2, b3)
    cf2 = C.fast_speed(rho, p, g1, b2, b3, b1)
    cf3 = C.fast_speed(rho, p, g1, b3, b1, b2)
    dv1 = np.abs(v1) + cf1
    dv2 = np.abs(v2) + cf2 / sn
    dv3 = np.abs(v3) + cf3 / sn
    dt1 = g['dx1'][None, None, None, :] / dv1
    dt2 = g['dx2'] / dv2
    dt3 = g['dx3'] / dv3
    return np.minimum(np.minimum(dt1, dt2), dt3), dt1, dt2, dt3, cf1, cf2, cf3


def report(name, path, fh, ntop=10):
    raw, d, g, geo, T, p, g1 = load(path)
    dtc, dt1, dt2, dt3, cf1, cf2, cf3 = cell_dt(d, g, T, p, g1)
    rho = d['dens']
    vabs = np.sqrt(d['velx']**2 + d['vely']**2 + d['velz']**2)
    bsq = d['bcc1']**2 + d['bcc2']**2 + d['bcc3']**2
    beta = np.where(bsq > 0, 2 * p / bsq, np.nan)
    rc = g['rc']
    lat = np.degrees(geo['lat'])
    lon = np.degrees(geo['lon'])
    flat = dtc.ravel()
    order = np.argsort(flat)[:ntop]
    idx = np.unravel_index(order, dtc.shape)
    fh.write('\n#### %s  (t = %.6e, cycle %d, dt_cfl_min = %.4f s)\n\n'
             % (name, raw['time'], raw['cycle'], C.CFL * flat[order[0]]))
    # which direction sets it, globally
    mins = [C.CFL * dt1.min(), C.CFL * dt2.min(), C.CFL * dt3.min()]
    fh.write('global per-direction limits: dt1(radial) %.3f  dt2(xi) %.3f  '
             'dt3(eta) %.3f s\n\n' % tuple(mins))
    fh.write('| # | cfl*dt [s] | dir | m,k,j,i | r/Rp | lat | lon_ss | p [bar] | T [K] '
             '| rho | on dfloor | beta | \\|v\\| | cf1 | dx_min |\n')
    fh.write('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n')
    for n in range(ntop):
        m, k, j, i = [int(x[n]) for x in idx]
        dirs = np.argmin([dt1[m, k, j, i], dt2[m, k, j, i], dt3[m, k, j, i]])
        dxm = [g['dx1'][i], g['dx2'][m, k, j, i], g['dx3'][m, k, j, i]][dirs]
        fh.write('| %d | %.4f | %s | %d,%d,%d,%d | %.3f | %+.1f | %+.1f | %.3e | %.0f '
                 '| %.3e | %s | %.2e | %.2e | %.2e | %.2e |\n'
                 % (n + 1, C.CFL * dtc[m, k, j, i], 'r xi eta'.split()[dirs], m, k, j, i,
                    rc[i] / 9.44e9, lat[m, k, j], lon[m, k, j], p[m, k, j, i] / 1e6,
                    T[m, k, j, i], rho[m, k, j, i],
                    'yes' if rho[m, k, j, i] <= 1.0001 * C.DFLOOR else 'no',
                    beta[m, k, j, i], vabs[m, k, j, i], cf1[m, k, j, i], dxm))
    # floor census
    onfl = rho <= 1.0001 * C.DFLOOR
    onpf = p <= 1.0001 * C.PFLOOR
    nz = onfl.sum(axis=(0, 1, 2))
    i0 = int(np.argmax(nz > 0)) if nz.any() else -1
    fh.write('\nfloored cells in this dump: rho<=dfloor %d (%.3f%%), p<=pfloor %d; '
             'lowest radial index with rho on dfloor i = %d (r/Rp %.3f)\n'
             % (onfl.sum(), 100 * onfl.mean(), onpf.sum(), i0,
                rc[i0] / 9.44e9 if i0 >= 0 else np.nan))
    if onfl.any():
        la = np.broadcast_to(lat[..., None], rho.shape)[onfl]
        lo = np.broadcast_to(lon[..., None], rho.shape)[onfl]
        day = np.abs(lo) < 90
        fh.write('  their |lat| mean %.1f deg, day-side fraction %.3f, '
                 'lon range %.0f..%.0f\n' % (np.abs(la).mean(), day.mean(),
                                             lo.min(), lo.max()))
    return dict(name=name, raw=raw, d=d, g=g, geo=geo, T=T, p=p, dtc=dtc,
                onfl=onfl, rho=rho)


if __name__ == '__main__':
    fh = open(C.OUT + 'cfl_tables.md', 'w')
    fh.write('# Per-cell CFL reconstruction (cfl.py)\n')
    fh.write('\nReproduces `mhd_newdt.cpp` exactly: dt_i = dx_i/(|v_i| + cf_i), with\n'
             'cf_2,3 divided by sin_cell, dx from the polynomial radial stretch and the\n'
             'gnomonic arc lengths, Gamma_1 and p from the run\'s own EOS table.\n')
    res = {}
    for name, dd in C.ARMS1:
        f = sorted(glob.glob(dd + '/bin/*.bin'))[-1]
        res[name] = report(name, f, fh)
    fh.close()
    print(open(C.OUT + 'cfl_tables.md').read())
