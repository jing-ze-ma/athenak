"""rstinfo.py <rst>: time, dt, ncycle of a restart file (restart.cpp STEP 1 layout)."""
import struct
import sys

d = open(sys.argv[1], 'rb').read(400000)
i = d.find(b'<par_end>')
i = d.find(b'\n', i) + 1
o = i + 4 + 4 + 9*8 + 2*19*4   # RegionIndcs = 19 ints
t, dt = struct.unpack('<dd', d[o:o+16])
nc = struct.unpack('<i', d[o+16:o+20])[0]
print('%.17e %.17e %d' % (t, dt, nc))
