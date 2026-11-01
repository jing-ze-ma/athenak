# am_vector spin reservoir from a hydro dhj restart: per-MeshBlock record = m slabs of (116,20,20) (m from the file size),
# the spin is the last slab of each record (restart.cpp kAmSpinRstMagic); 24 MeshBlocks.
import numpy as np, os, sys
f = sys.argv[1]; n = 116*20*20; sz = os.path.getsize(f); m = sz // (24*n*8); R = m*n*8
tot = 0.0; nz = 0; mx = 0.0
for k in range(24):
    a = np.fromfile(f, dtype='<f8', count=n, offset=sz - n*8 - k*R)
    nz += np.count_nonzero(a); tot += a.sum(); mx = max(mx, np.abs(a).max())
print(f"{f}: slabs/record {m}, spin nonzero cells {nz} of {24*112*16*16} active, sum {tot:.4e}, max|s| {mx:.4e}")
