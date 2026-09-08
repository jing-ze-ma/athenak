"""The gate an open lower boundary has to pass: does the star radiate what it is fed?

With problem/inner_bc = open the energy enters as the ENTROPY of the inflow at the lowest
active layer, so the luminosity is an OUTPUT.  This compares the net longwave flux the
two-stream carries at the top of a column against L/(4 pi r^2), before and after
relaxation, using the solver's own column dumps (problem/ck_dump_file and ck_dump_file2).

usage: flux_gate.py L_erg_per_s COLUMN_DUMP [COLUMN_DUMP_LATE]
"""
import sys
import numpy as np

L = float(sys.argv[1])
SIG = 5.670374419e-5


def report(fn, label):
    c = np.loadtxt(fn)
    i, rf, pbar, T, Flw, Qsw, g1, gad, tau, w = c.T
    k = int(np.argmax(rf))
    floc = L/(4.0*np.pi*rf[k]**2)
    teff = (Flw[k]/SIG)**0.25 if Flw[k] > 0 else float('nan')
    print('  %-6s top face r = %.4e cm : F_lw = %.4e, L/(4 pi r^2) = %.4e'
          % (label, rf[k], Flw[k], floc))
    print('         F_lw / F_local = %.3f   -> effective temperature %.0f K'
          % (Flw[k]/floc, teff))
    # where the two-stream still owns a share, and what it carries there
    m = w < 0.999
    if m.sum():
        fl = Flw[m]/(L/(4.0*np.pi*rf[m]**2))
        print('         over the %d faces the RT owns: F/F_local = %.3f .. %.3f'
              % (m.sum(), fl.min(), fl.max()))
    return Flw[k]/floc


a = report(sys.argv[2], 'start')
if len(sys.argv) > 3:
    b = report(sys.argv[3], 'late')
    print()
    print('  the gate: %.3f -> %.3f of the required flux%s'
          % (a, b, '  (moving toward 1)' if abs(b-1) < abs(a-1) else '  (NOT improving)'))
