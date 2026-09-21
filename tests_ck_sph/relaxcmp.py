"""Gate (e): one radial column of the relaxation runs, switch off vs on.

The binary output of this problem carries (dens, velx, vely, velz, eint) only -- under the
tabulated EOS there is no temperature in it -- so what is compared here is the SPECIFIC
internal energy e/rho, which is a monotone proxy for T at fixed composition, and the
density.  The radius and the t = 0 pressure of each cell are taken from the gate-(c)
column dump, which is on the same (stretched) radial grid.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "vis", "python"))
import bin_convert                                            # noqa: E402

KK, JJ, MB = 3, 5, 0


def column(run):
    files = sorted(os.listdir(os.path.join(run, "bin")))
    d = bin_convert.read_binary(os.path.join(run, "bin", files[-1]))
    return (d["time"], np.array(d["mb_data"]["dens"])[MB, KK, JJ, :],
            np.array(d["mb_data"]["eint"])[MB, KK, JJ, :], files[-1])


if __name__ == "__main__":
    ref = np.loadtxt("c_deep_off/col.txt")
    rf, p0 = ref[:, 1], ref[:, 2]
    rc = 0.5*(rf[:-1] + rf[1:])                       # cell centres, same grid
    t0, d0, e0, f0 = column("e_relax_off")
    t1, d1, e1, f1 = column("e_relax_on")
    print(f"off: {f0}  t = {t0:.5e} s      on: {f1}  t = {t1:.5e} s")
    n = d0.size
    rc, p0 = rc[:n], p0[:n]
    s0, s1 = e0/d0, e1/d1
    rel = s1/s0 - 1.0
    print("  i    r[cm]     p(t=0)[bar]   e/rho off    e/rho on     rel     "
          "rho off      rho on      rel")
    for i in range(0, n, max(1, n//26)):
        print(f"{i:4d} {rc[i]:.4e} {p0[i]:.4e} {s0[i]:.5e} {s1[i]:.5e} "
              f"{rel[i]:+.4f}  {d0[i]:.4e} {d1[i]:.4e} {d1[i]/d0[i]-1:+.4f}")
    j = int(np.argmax(np.abs(rel)))
    print(f"max |d(e/rho)|/(e/rho) = {abs(rel[j]):.4f} at i = {j}, "
          f"r = {rc[j]:.4e}, p(t=0) = {p0[j]:.3e} bar")
    # near the photosphere band, 1e-3 .. 1 bar
    band = (p0 > 1.0e-3) & (p0 < 1.0)
    if band.any():
        print(f"photospheric band 1e-3..1 bar: mean rel {rel[band].mean():+.4f}, "
              f"min {rel[band].min():+.4f}, max {rel[band].max():+.4f}")
