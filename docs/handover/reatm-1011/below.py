"""below.py RPH runA runB: per common dump, the outermost radial cell index (and r/R_ph) up to which runA and runB
are BITWISE identical in every hydro and m1 variable (all theta/phi), and the max rel diff in r < R_ph and r < 0.99 R_ph."""
import sys, glob, os, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import bin_convert as bc
R = float(sys.argv[1]); a, b = sys.argv[2], sys.argv[3]
for kind in ('hydro_w', 'm1'):
    for f in sorted(glob.glob(a + '/bin/agcar3d.%s.*.bin' % kind)):
        g = f.replace(a, b, 1)
        if not os.path.exists(g): continue
        ha, hb = bc.read_binary_as_athdf(f), bc.read_binary_as_athdf(g); c = int(ha['NumCycles'])
        r = np.asarray(ha['x1v'], float)/R; first = len(r); out = []
        for v in [k for k in ha if k not in ('x1f','x1v','x2f','x2v','x3f','x3v','Time','NumCycles','MaxLevel')]:
            x, y = np.asarray(ha[v], float), np.asarray(hb[v], float)
            dd = (np.abs(y - x)).max(axis=(0, 1)); sc = np.abs(x).max(axis=(0, 1)) + 1e-300
            nz = np.where(dd > 0)[0]
            if len(nz): first = min(first, nz[0])
            rel = dd/sc
            out.append('%s %.1e/%.1e' % (v, rel[r < 1.0].max(), rel[r < 0.99].max()))
        print('c%-4d %-7s first differing cell i=%d (r/R_ph %.4f)  max rel diff r<1 / r<0.99: %s' % (
            c, kind, first, r[first] if first < len(r) else -1, '  '.join(out)))
