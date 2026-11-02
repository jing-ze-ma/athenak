"""
Advected electron fraction <mhd>/xe_tracer (xetracer-1101) on a uniform MHD box with the
general EOS table (src/pgen/xe_tracer_test.cpp, inputs/tests/xe_tracer_box.athinput).

1. Relaxation: z0 = 10 x_eq (and 0.1 x_eq) at rest; after N cycles z must equal the
   backward-Euler recursion n <- 2c/(1 + sqrt(1 + 4ac)), c = n + a n_eq^2,
   a = alpha_Na(T) dt (Verner & Ferland 1996, recomputed here), to round-off.  Measured
   on xetracer-1101: rel. diff 1.8e-15 / 1.1e-16.
2. Advection (xe_alpha_fac = 0): a Gaussian bump in z carried once across the periodic
   box at nx1 = 64, 128.  L1(z - z0)/L1(bump) measured 0.410 / 0.158; sum z conserved to
   1e-14; no undershoot below the background.
3. Restart: a file written with the tracer continues bitwise (final restart cmp); a file
   written without it (xe_tracer = false) restarted with mhd/xe_tracer=true starts from
   local Saha (z = x_eq exactly at rest).
"""

import filecmp
import os
import subprocess

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(REPO, "inputs", "tests", "xe_tracer_box.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_xe_tracer")

RELAX_TOL = 1.0e-12
ADV_L1_MAX = {64: 0.45, 128: 0.18}
ADV_RATIO_MIN = 2.3       # measured 2.60


def build():
    exe = os.path.join(BUILD, "src", "athena")
    if os.path.exists(exe):
        return exe
    subprocess.run(["cmake", "-S", REPO, "-B", BUILD, "-D", "PROBLEM=xe_tracer_test",
                    "-D", "CMAKE_BUILD_TYPE=Release"],
                   check=True, capture_output=True, text=True)
    subprocess.run(["make", "-C", BUILD, "-j", "16"],
                   check=True, capture_output=True, text=True)
    return exe


def alpha_na(t):
    a, b, t0, t1 = 5.641e-12, 0.1749, 3.077e2, 2.617e6
    s0, s1 = np.sqrt(t/t0), np.sqrt(t/t1)
    return a/(s0*(1 + s0)**(1 - b)*(1 + s1)**(1 + b))


def run(exe, tag, args, restart=None):
    d = os.path.join(BUILD, "runs", tag)
    os.makedirs(d, exist_ok=True)
    cmd = [exe, "-r", restart] if restart else [exe, "-i", INPUT]
    proc = subprocess.run(cmd + ["-d", d] + args, cwd=d, capture_output=True, text=True)
    if proc.returncode != 0 or "FATAL" in proc.stdout:
        pytest.fail(f"{tag} failed:\n{proc.stdout[-3000:]}\n{proc.stderr[-2000:]}")
    rows = [[float(v) for v in ln.split()[1:]] for ln in proc.stdout.splitlines()
            if ln.startswith("XE ")]
    return d, proc.stdout, np.array(rows)


def test_xe_tracer_cpu():
    exe = build()
    # 1. relaxation from above and from below
    for zf in (10.0, 0.1):
        _, _, r = run(exe, f"relax{zf}", [f"problem/z_fac={zf}", "time/nlim=10"])
        z, xeq, nt, t, dt, nc = r[:, 2], r[0, 3], r[0, 4], r[0, 5], r[0, 6], int(r[0, 7])
        a, neq, n = alpha_na(t)*dt, xeq*nt, zf*xeq*nt
        for _ in range(nc):
            c = n + a*neq**2
            n = 2*c/(1 + np.sqrt(1 + 4*a*c))
        err = np.max(np.abs(z - n/nt))/(n/nt)
        print(f"relax z0/x_eq={zf}: z/x_eq={z[0]/xeq:.15g} BE={n/nt/xeq:.15g} "
              f"err {err:.2e}")
        assert err < RELAX_TOL, f"relaxation z0/x_eq={zf}: rel diff {err:.3e}"
    # 2. advection
    l1 = {}
    for nx in (64, 128):
        _, _, r = run(exe, f"adv{nx}", [
            "problem/z_fac=1", "problem/z_amp=5", "problem/z_width=0.05",
            "problem/vx=1.0e6", "mhd/xe_alpha_fac=0.0", "time/nlim=-1",
            "time/tlim=1000.0", f"mesh/nx1={nx}", f"meshblock/nx1={nx}"])
        x, z, xeq = r[:, 1], r[:, 2], r[:, 3]
        z0 = xeq*(1 + 5*np.exp(-((x - 0.5e9)/0.05e9)**2))
        l1[nx] = np.abs(z - z0).sum()/np.abs(z0 - xeq).sum()
        cons = abs(z.sum() - z0.sum())/z0.sum()
        print(f"advect nx={nx}: L1 {l1[nx]:.4e} sum z {cons:.2e} "
              f"min {z.min()/xeq[0]:.6f}")
        assert l1[nx] < ADV_L1_MAX[nx], f"advection L1 at nx={nx}: {l1[nx]:.3e}"
        assert cons < 1.0e-12, f"scalar not conserved at nx={nx}: {cons:.3e}"
        assert z.min() >= xeq[0]*(1 - 1.0e-12), f"undershoot at nx={nx}"
    assert l1[64]/l1[128] > ADV_RATIO_MIN, f"advection order: {l1[64]/l1[128]:.2f}"
    # 3. restarts
    k = ["problem/z_fac=3", "problem/z_amp=5", "problem/z_width=0.05",
         "problem/vx=1.0e6", "mhd/xe_alpha_fac=0.1"]
    ds, _, _ = run(exe, "rst_straight", k + ["time/nlim=20"])
    dc, _, _ = run(exe, "rst_cont", ["time/nlim=20"],
                   restart=os.path.join(ds, "rst", "xebox.00001.rst"))
    assert filecmp.cmp(os.path.join(ds, "rst", "xebox.00002.rst"),
                       os.path.join(dc, "rst", "xebox.00002.rst"), shallow=False), \
        "restart with the tracer is not a bitwise continuation"
    dn, _, _ = run(exe, "rst_notr", k + ["mhd/xe_tracer=false", "time/nlim=10"])
    _, out, r = run(exe, "rst_on", ["mhd/xe_tracer=true", "time/nlim=11"],
                    restart=os.path.join(dn, "rst", "xebox.00001.rst"))
    assert "starts from local Saha" in out
    assert np.max(np.abs(r[:, 2]/r[:, 3] - 1)) < 1.0e-12, "switch-on z != x_eq"
    # 4. density floor: a floored cell keeps its y = n_e/n_H (the floor mass inherits it)
    _, out, r = run(exe, "floor", ["problem/z_fac=1", "mhd/xe_alpha_fac=0.0",
                                   "mhd/dfloor=2.0e-8", "time/nlim=2"])
    y0 = float([ln for ln in out.splitlines() if ln.startswith("XEINIT")][0].split()[-1])
    yerr = np.max(np.abs(r[:, 9]/y0 - 1))
    print(f"floor: y/y0 - 1 max {yerr:.2e}")
    assert yerr < 1.0e-12, f"floored cell changed y: {yerr:.3e}"
