"""cmphst.py A.hst B.hst [t0]: rows of B (data rows, t >= t0) that are also time rows of A: count verbatim-equal vs differing.
Use: BASE window (B) vs the production run's hst over the window (production read-only): expected all verbatim (GATE-BASE extended)."""
import sys
def rows(f):
    d = {}
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        d[l.split()[0]] = l.rstrip('\n')          # last occurrence wins (restart overlaps)
    return d
A, B = rows(sys.argv[1]), rows(sys.argv[2]); t0 = float(sys.argv[3]) if len(sys.argv) > 3 else -1e300
ks = [k for k in B if float(k) >= t0 and k in A]
bad = [k for k in ks if A[k] != B[k]]
print('%s vs %s: %d common rows (t >= %g), %d verbatim, %d differ%s' % (sys.argv[1], sys.argv[2], len(ks), t0, len(ks) - len(bad),
      len(bad), ('; first differing t = ' + min(bad, key=float)) if bad else ''))
print('HSTCMP ' + ('BITWISE' if ks and not bad else ('NO COMMON ROWS' if not ks else 'DIFF')))
