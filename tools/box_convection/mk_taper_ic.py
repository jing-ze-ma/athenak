"""Re-express a relaxed column in the TAPERED EOS at fixed (rho, T).

Turning the thin-region radiation taper on for a column that was relaxed without
it needs this first.  The input lines the taper adds, for the 15 Msun B-star box:

    <hydro>
    eos_radiation  = true
    eos_rad_rho_hi = 1.0e-8     # w = 1 (full LTE aT^4) at and below this density
    eos_rad_rho_lo = 1.0e-9     # w = 0 (the two-stream owns the radiation) above
    <problem>
    rt_rad_force   = true       # the momentum source that goes with the taper

and the EOS-table dump this script reads comes from a zero-cycle run with
<hydro>/eos_table_dump set.


The taper changes what e MEANS: e = e_gas(rho,T) + w(rho) a T^4 instead of
e_gas + a T^4.  Feeding the untapered column's e straight to the tapered EOS
therefore preserves e and throws all of it into the gas, which at the top --
where a T^4 is 17x the gas thermal energy -- multiplies T by ~7 and starts the
run with a 40x F_bot inward flux.  Preserving (rho, T) instead keeps the
radiation field and the temperature structure and changes only the energy
BOOK-KEEPING, which is exactly what the taper is.

Usage: python3 mk_taper_ic.py <eos_table_dump> <in_ic.txt> <out_ic.txt> rho_lo rho_hi
"""
import sys
import numpy as np

A_RAD = 7.5657332503e-15


def read_dump(fn):
    with open(fn) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    le = d[:, 0].reshape(ny, nx)
    lp = d[:, 1].reshape(ny, nx)
    x = xmin + dx*np.arange(nx)
    y = ymin + dy*np.arange(ny)
    return x, y, le, lp


def bilin(x, y, f, xq, yq):
    ix = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    iy = np.clip(np.searchsorted(y, yq) - 1, 0, len(y) - 2)
    tx = (xq - x[ix])/(x[ix+1] - x[ix])
    ty = (yq - y[iy])/(y[iy+1] - y[iy])
    return ((1-tx)*(1-ty)*f[iy, ix] + tx*(1-ty)*f[iy, ix+1]
            + (1-tx)*ty*f[iy+1, ix] + tx*ty*f[iy+1, ix+1])


def weight(lrho, llo, lhi):
    s = np.clip((lrho - llo)/(lhi - llo), 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def main():
    dump, fin, fout, rlo, rhi = sys.argv[1:6]
    rlo, rhi = float(rlo), float(rhi)
    x, y, le, lp = read_dump(dump)
    a = np.loadtxt(fin, comments='#')
    z, rho, e = a[:, 0], a[:, 1], a[:, 2]
    lr = np.log10(rho)
    # invert e = rho*10^E(lr,lT) + a T^4 for T, with the UNTAPERED EOS
    lt = np.full_like(lr, 4.7)
    for _ in range(80):
        eg = rho*10.0**bilin(x, y, le, lr, lt)
        f = eg + A_RAD*10.0**(4*lt) - e
        d = 1.0e-4
        eg2 = rho*10.0**bilin(x, y, le, lr, lt + d)
        df = (eg2 + A_RAD*10.0**(4*(lt + d)) - e - f)/d
        lt = np.clip(lt - f/df, y[0], y[-1])
    T = 10.0**lt
    eg = rho*10.0**bilin(x, y, le, lr, lt)
    w = weight(lr, np.log10(rlo), np.log10(rhi))
    enew = eg + w*A_RAD*T**4
    rel = np.abs(eg + A_RAD*T**4 - e)/e
    print('T inversion max rel residual %.3e' % rel.max())
    print('%12s %12s %12s %8s %10s %10s' % ('z', 'rho', 'T', 'w', 'e_old',
                                            'e_new'))
    for i in range(0, len(z), max(1, len(z)//14)):
        print('%12.4e %12.4e %12.4e %8.4f %10.4e %10.4e'
              % (z[i], rho[i], T[i], w[i], e[i], enew[i]))
    hdr = ('# bstar_fecz relaxed column, re-expressed in the TAPERED EOS\n'
           '# (<hydro>/eos_rad_rho_hi = %.3g, eos_rad_rho_lo = %.3g) at FIXED\n'
           '# (rho, T): e = e_gas(rho,T) + w(rho) a T^4.  Source: %s\n'
           '# z[cm]  rho[g/cm^3]  eint[erg/cm^3]\n' % (rhi, rlo, fin))
    with open(fout, 'w') as fh:
        fh.write(hdr)
        for zz, rr, uu in zip(z, rho, enew):
            fh.write('%.10e %.10e %.10e\n' % (zz, rr, uu))
    print('wrote', fout, len(z), 'nodes')


main()
