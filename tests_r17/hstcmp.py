import sys
import numpy as np
TURN = 4705.0
for arm in sys.argv[1:]:
    d = np.loadtxt(arm + '/he4.hydro.hst')
    out = []
    for f in (0.25, 0.5, 1.0, 1.5):
        k = int(np.argmin(np.abs(d[:, 0] - f*TURN)))
        if abs(d[k, 0] - f*TURN) > 0.1*TURN:
            out.append("%.2f: --" % f)
            continue
        out.append("%.2f: mom %.3e KEr %.3e" % (f, abs(d[k, 3]), d[k, 7]))
    print("%-22s %s" % (arm, " | ".join(out)))
