# rstcyc.py <rst> <dt>: print ncycle of an AthenaK restart (header: ... Real time, Real dt, int ncycle after the input text)
import sys, struct
b = open(sys.argv[1], 'rb').read(4 << 20); dt = float(sys.argv[2])
i = b.find(b'<par_end>'); i = 0 if i < 0 else i
for o in range(i, len(b) - 20):
    t, d = struct.unpack_from('<dd', b, o)
    if abs(d - dt) < 1e-3 * dt and 0 <= t < 1e9:
        n = struct.unpack_from('<i', b, o + 16)[0]
        if 0 <= n < 10**7 and abs(t - n * dt) < 1e-6 * max(t, dt): print(n); sys.exit(0)
print(-1)
