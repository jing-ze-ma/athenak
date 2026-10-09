#!/usr/bin/env python3
"""Build a TOPS mixture string (NUMBER fractions) for a GS98 metal mixture at
given X (H mass fraction) and Z (metal mass fraction).  Metal number fractions
are the Grevesse-Sauval (1998) 'standard mix' as published by the LANL
astrophysical opacity page, https://aphysics2.lanl.gov/static/opacdocs/grevsau1.html
"""
import sys

# Z, symbol, number fraction within the metals (GS98, LANL grevsau1 page)
GS98 = [
    (6, 'c', 2.45187e-01), (7, 'n', 6.15882e-02), (8, 'o', 5.00608e-01),
    (9, 'f', 2.68842e-05), (10, 'ne', 8.90220e-02), (11, 'na', 1.58306e-03),
    (12, 'mg', 2.81512e-02), (13, 'al', 2.18523e-03), (14, 'si', 2.62723e-02),
    (15, 'p', 2.08688e-04), (16, 's', 1.58306e-02), (17, 'cl', 2.34152e-04),
    (18, 'ar', 1.85993e-03), (19, 'k', 9.76107e-05), (20, 'ca', 1.69628e-03),
    (21, 'sc', 1.09521e-06), (22, 'ti', 7.75349e-05), (23, 'v', 7.40453e-06),
    (24, 'cr', 3.46336e-04), (25, 'mn', 1.81760e-04), (26, 'fe', 2.34152e-02),
    (27, 'co', 6.15882e-05), (28, 'ni', 1.31673e-03), (29, 'cu', 1.20087e-05),
    (30, 'zn', 2.94780e-05)]
# standard atomic weights
A = {'c': 12.011, 'n': 14.007, 'o': 15.999, 'f': 18.998, 'ne': 20.180,
     'na': 22.990, 'mg': 24.305, 'al': 26.982, 'si': 28.085, 'p': 30.974,
     's': 32.06, 'cl': 35.45, 'ar': 39.948, 'k': 39.098, 'ca': 40.078,
     'sc': 44.956, 'ti': 47.867, 'v': 50.942, 'cr': 51.996, 'mn': 54.938,
     'fe': 55.845, 'co': 58.933, 'ni': 58.693, 'cu': 63.546, 'zn': 65.38,
     'h': 1.008, 'he': 4.0026}


def mixture(X, Z, drop_below=1e-7):
    """Return (mixture string, dict of realised mass fractions)."""
    Y = 1.0 - X - Z
    assert Y > -1e-12, 'X+Z>1'
    # mean metal weight, so that a metal number density can be converted
    mbar_Z = sum(f*A[s] for _, s, f in GS98)
    # per unit total mass: n_i propto (mass fraction)/(weight)
    n = {}
    if X > 0.0:
        n['h'] = X/A['h']
    if Y > 0.0:
        n['he'] = Y/A['he']
    for _, s, f in GS98:
        n[s] = Z*f/mbar_Z          # f/mbar_Z = number per unit metal mass
    tot = sum(n.values())
    n = {k: v/tot for k, v in n.items()}
    n = {k: v for k, v in n.items() if v >= drop_below}
    tot = sum(n.values())
    n = {k: v/tot for k, v in n.items()}
    mt = sum(v*A[k] for k, v in n.items())
    mass = {k: v*A[k]/mt for k, v in n.items()}
    s = ' '.join('%.5e %s' % (n[k], k) for k in n)
    return s, mass


if __name__ == '__main__':
    X, Z = float(sys.argv[1]), float(sys.argv[2])
    s, mass = mixture(X, Z)
    print(len(s), 'chars')
    print(s)
    print('X=%.5f Y=%.5f Z=%.5f' % (mass.get('h', 0.0), mass.get('he', 0.0),
                                    sum(v for k, v in mass.items()
                                        if k not in ('h', 'he'))))
    print('Fe mass frac %.3e' % mass.get('fe', 0.0))
