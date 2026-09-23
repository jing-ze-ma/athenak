# print time and ncycle of an AthenaK restart: after <par_end>, nmb,rootlev,RegionSize,
# RegionIndcs x2, time, dt, ncycle; find time/dt as the two doubles before ncycle by
# scanning for a plausible (time>0, 0<dt<time) pair followed by an int.
import sys, struct
b = open(sys.argv[1], 'rb').read(200000)
p = b.index(b'<par_end>') + len(b'<par_end>')
while b[p:p+1] in (b'\n', b' '): p += 1
for o in range(p + 8, p + 400):
    t, dt = struct.unpack('<dd', b[o:o+16]); n = struct.unpack('<i', b[o+16:o+20])[0]
    if 1e3 < t < 1e12 and 1e-3 < dt < 1e4 and 1000 < n < 10**8:
        print('time %.6e dt %.4e ncycle %d offset %d' % (t, dt, n, o - p)); break
