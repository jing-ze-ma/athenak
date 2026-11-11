"""raycheap: bitwise compare of the last m1 / m1_lat bins of two run dirs (MeshBlocks sorted by geometry, so different
rank layouts compare).  python3 cmpbin.py RUN_A RUN_B   (needs vis/python/bin_convert.py of this repo on sys.path)"""
import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'vis', 'python'))
import bin_convert as bc
def ld(d, v):
    f = sorted(glob.glob(f"{d}/bin/*.{v}.*.bin"))[-1]; x = bc.read_binary(f)
    o = sorted(range(len(x["mb_geometry"])), key=lambda m: tuple(np.round(x["mb_geometry"][m], 12)))
    return x["cycle"], {k: np.concatenate([np.asarray(x["mb_data"][k][m]).ravel() for m in o]) for k in x["var_names"]}
for v in ("m1", "m1_lat"):
    ca, A = ld(sys.argv[1], v); cb, B = ld(sys.argv[2], v)
    nb = sum(int(np.array_equal(A[k], B[k])) for k in A)
    dm = max(float(np.max(np.abs(A[k] - B[k]))) for k in A)
    print(f"{v}: cycles {ca}/{cb} bitwise {nb}/{len(A)} vars, max |diff| {dm:.3e}")
