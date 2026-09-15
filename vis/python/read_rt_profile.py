"""Reader for the box_convection horizontally averaged x1 profile dump.

The dump is written by src/pgen/box_convection.cpp when problem/rt_profile_dt > 0
(file name: problem/rt_profile_file, default "rt_profile.bin").  It is a stream of
records, appended once per cycle at the firing cadence, each

    [float64 time][int32 nx1][int32 nvar][float64 x1v[nx1]][float64 data[nvar][nx1]]

little-endian, no padding, no global header.  Every row of `data` is the mean over
the WHOLE horizontal plane (all x2, all x3, all MeshBlocks, all MPI ranks) at that
x1 index, taken from the primitives w0 (plus u0(IEN) and Hydro::wtemp) at the start
of the cycle, in code units.  The nvar = 8 rows are named in VARS below.
"""

import numpy as np

VARS = ["rho", "v1", "rhov1", "v1sq", "vhsq", "T", "e", "fenth"]


def read_rt_profile(fname, last_only=False):
    """Read rt_profile.bin.

    Returns a dict with 'time' (nrec,), 'x1v' (nx1,) and one (nrec, nx1) array per
    entry of VARS.  With last_only=True the record axis is dropped and only the
    final record is returned.
    """
    times, x1v, recs = [], None, []
    with open(fname, "rb") as f:
        while True:
            head = f.read(16)
            if len(head) < 16:
                break
            t = np.frombuffer(head[:8], "<f8")[0]
            nx1, nvar = np.frombuffer(head[8:], "<i4")
            nx1, nvar = int(nx1), int(nvar)
            x = np.frombuffer(f.read(8*nx1), "<f8")
            d = np.frombuffer(f.read(8*nvar*nx1), "<f8").reshape(nvar, nx1)
            if x1v is None:
                x1v = x.copy()
            times.append(t)
            recs.append(d.copy())
    if not recs:
        raise IOError("no records in " + fname)
    arr = np.array(recs)                       # (nrec, nvar, nx1)
    out = {"time": np.array(times), "x1v": x1v}
    for n, name in enumerate(VARS[:arr.shape[1]]):
        out[name] = arr[-1, n] if last_only else arr[:, n]
    if last_only:
        out["time"] = out["time"][-1]
    return out
