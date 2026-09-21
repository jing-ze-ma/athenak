"""Gate (c): the chord geometry of the pseudo-spherical beam, checked two ways.

1.  The closed-form chord the kernel uses,
        ds_j = sqrt(r_{j+1}^2 - b^2) - sqrt(max(r_j^2 - b^2, 0)),  b = r sin(theta0),
    with weight 2 on the shells below the target when mu0 < 0 (near + far leg), against a
    DIRECT numerical integration of the straight ray x(t) = P + t s_hat: march t in small
    steps, bin each step into the shell that |x(t)| falls in, and compare the binned path
    lengths.  Gate 1e-6 relative.

2.  tau_ray for an isothermal exponential atmosphere of constant opacity against the
    CHAPMAN function Ch(X, theta0), X = r/H, computed by its own quadrature.  This is the
    textbook closed form the pseudo-spherical approximation is supposed to reproduce.
"""
import numpy as np

# production radial grid
R0, R1, NR = 9.44e9, 2.0556e10, 128
RF = np.linspace(R0, R1, NR + 1)


def chords(r, mu0, rf, rcut):
    """Path lengths per shell for a ray to the target at radius r.  Returns None if dark."""
    b = r*np.sqrt(max(0.0, 1.0 - mu0*mu0))
    n = rf.size - 1
    ds = np.zeros(n)
    # index of the target's own shell (its lower face is rf[f] == r when r is a face)
    f = int(np.searchsorted(rf, r, side="right") - 1)
    f = min(max(f, 0), n - 1)
    if mu0 >= 0.0:
        jlo = f
    else:
        if b <= rcut:
            return None
        jlo = int(np.searchsorted(rf, b, side="right") - 1)
        jlo = max(jlo, 0)
    prev = np.sqrt(max(rf[jlo]**2 - b*b, 0.0))
    if mu0 >= 0.0:
        prev = np.sqrt(max(r*r - b*b, 0.0))          # start at the target itself
    for j in range(jlo, n):
        cur = np.sqrt(max(rf[j+1]**2 - b*b, 0.0))
        w = 2.0 if j < f else 1.0
        ds[j] = w*(cur - prev)
        prev = cur
    if mu0 < 0.0 and f < n:
        # the target's own shell: near leg only from r down to its lower face
        cur = np.sqrt(max(rf[f+1]**2 - b*b, 0.0))
        lo = np.sqrt(max(rf[f]**2 - b*b, 0.0))
        ds[f] = (cur - lo) + (r*abs(mu0) - lo)
    return ds


def chords_numeric(r, mu0, rf):
    """Independent reference: for every shell face solve |P + t s_hat| = r_f for t with a
    polynomial root finder (np.roots on t^2 + 2 t r mu0 + r^2 - r_f^2), collect all
    crossings along t > 0, sort them, and attribute each interval between consecutive
    crossings to the shell its midpoint falls in.  No sqrt identity of the kernel's is
    reused, and no marching error."""
    n = rf.size - 1
    ts = [0.0]
    for rfj in rf:
        c = r*r - rfj*rfj
        rt = np.roots([1.0, 2.0*r*mu0, c])
        for z in rt:
            if abs(z.imag) < 1.0e-9*max(1.0, abs(z.real)) and z.real > 1.0e-9:
                ts.append(float(z.real))
    ts = np.unique(np.array(sorted(ts)))
    ds = np.zeros(n)
    for a, b in zip(ts[:-1], ts[1:]):
        tm = 0.5*(a + b)
        rm = np.sqrt(r*r + 2.0*tm*r*mu0 + tm*tm)
        k = int(np.searchsorted(rf, rm, side="right") - 1)
        if 0 <= k < n:
            ds[k] += (b - a)
    return ds


def geometry_gate():
    print("=== (c1) chord formula vs direct ray integration ===")
    worst = 0.0
    for mu0 in [1.0, 0.9787, 0.4347, 0.1020, -0.1020, -0.1947, -0.5128]:
        for ir in [10, 40, 80, 120]:
            r = RF[ir]
            a = chords(r, mu0, RF, RF[0])
            if a is None:
                continue
            b = chords_numeric(r, mu0, RF)
            m = a > 1.0e-3*a.max()
            rel = np.max(np.abs(a[m] - b[m])/a[m])
            worst = max(worst, rel)
        print(f"  mu0 = {mu0:+.4f}: worst relative chord error over 4 target radii "
              f"= {worst:.3e}")
    print(f"  WORST OVERALL = {worst:.3e}")
    # mu0 = 1 must give exactly dr
    a = chords(RF[10], 1.0, RF, RF[0])
    dr = np.diff(RF)
    print(f"  mu0 = 1: max|ds_j - dr_j| over the shells above the target = "
          f"{np.max(np.abs(a[10:] - dr[10:])):.3e} cm (dr = {dr[0]:.4e})")


def chapman_gate():
    print("\n=== (c2) tau_ray vs the Chapman function, isothermal exponential atmosphere ===")
    print("    n(r) = n0 exp(-(r-r0)/H), constant kappa; tau_vert(r0) normalised to 1")
    for H in [1.24e9/8.0, 1.24e9]:
        r0 = RF[0]
        X = r0/H
        print(f"  H = {H:.3e} cm, X = r/H = {X:.1f}")
        for th in [60.0, 85.0, 90.0, 95.0]:
            mu0 = np.cos(np.radians(th))
            # discretised tau along the ray with the kernel's chords
            rc = 0.5*(RF[:-1] + RF[1:])
            kr = np.exp(-(rc - r0)/H)/H          # int kr dr from r0 to inf = 1
            rtar = r0 + 6.0*H     # high enough that a 95 deg ray's tangent point is
            ds = chords(rtar, mu0, RF, 0.0)   # still inside the grid, low enough to be in it
            tvert = np.exp(-(rtar - r0)/H)    # vertical depth above the target
            tau = np.nan if ds is None else float(np.sum(ds*kr))
            # Chapman by quadrature: tau(theta0)/tau_vert
            s = np.linspace(0.0, 200.0*H, 2000001)
            rr = np.sqrt(rtar*rtar + 2.0*s*rtar*mu0 + s*s)
            ch = float(np.trapezoid(np.exp(-(rr - r0)/H)/H, s))
            print(f"    theta0 = {th:5.1f} deg: Ch_kernel = {tau/tvert:9.4f}   "
                  f"Ch_exact = {ch/tvert:9.4f}   ratio = {tau/ch:8.5f}   "
                  f"secant 1/mu0 = {1.0/mu0 if mu0 > 0 else float('inf'):10.4f}")


if __name__ == "__main__":
    geometry_gate()
    chapman_gate()
