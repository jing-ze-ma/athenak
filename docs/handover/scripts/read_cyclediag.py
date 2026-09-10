#!/usr/bin/env python3
"""Reader for the per-cycle single-meshblock diagnostic dump ("cyclediag").

The dump is written by deep_hot_jupiter_rt's DhjCycleDiag when problem/diag_gid >= 0,
once per cycle, into cyclediag/<basename>.cyclediag.<cycle:08d>.dat.

Layout: a text header line, a text line of variable names, then raw little-endian
float64 arrays in that order.  Cell-centred variables are (n3, n2, n1) in C order,
rt_icut is (n3, n2), rad_w and rad_tauf are (n3, n2, n1+1) -- they live on x1 faces.
"""

import sys

import numpy as np

# variables that are not (n3, n2, n1)
FACE_VARS = ("rad_w", "rad_tauf")
COLUMN_VARS = ("rt_icut",)


def load(path):
    """Load a cyclediag file.  Returns a dict of the header fields plus one numpy
    array per variable name."""
    with open(path, "rb") as f:
        hdr = f.readline().decode().strip()
        names = f.readline().decode().strip().split()
        raw = f.read()

    if not hdr.startswith("cyclediag v1"):
        raise ValueError("%s: not a cyclediag v1 file (header %r)" % (path, hdr[:40]))

    out = {}
    for tok in hdr.split()[2:]:
        key, _, val = tok.partition("=")
        try:
            out[key] = int(val)
        except ValueError:
            out[key] = float(val)
    out["names"] = names

    n1, n2, n3 = out["n1"], out["n2"], out["n3"]
    nface = out["nface"]
    arr = np.frombuffer(raw, dtype="<f8")

    off = 0
    for name in names:
        if name in COLUMN_VARS:
            shape = (n3, n2)
        elif name in FACE_VARS:
            shape = (n3, n2, nface)
        else:
            shape = (n3, n2, n1)
        nel = int(np.prod(shape))
        if off + nel > arr.size:
            raise ValueError("%s: truncated at variable %r" % (path, name))
        out[name] = arr[off:off + nel].reshape(shape)
        off += nel
    if off != arr.size:
        raise ValueError("%s: %d trailing doubles after the last variable"
                         % (path, arr.size - off))
    return out


def main(argv):
    if len(argv) < 2:
        print("usage: read_cyclediag.py <file.dat> [k] [j]")
        return 1
    d = load(argv[1])
    k = int(argv[2]) if len(argv) > 2 else (d["ks"] + d["ke"]) // 2
    j = int(argv[3]) if len(argv) > 3 else (d["js"] + d["je"]) // 2
    is_, ie = d["is"], d["ie"]

    print("cycle %d  time %.10g  dt %.10g  gid %d  (k=%d, j=%d)"
          % (d["cycle"], d["time"], d["dt"], d["gid"], k, j))
    print("icut = %g" % d["rt_icut"][k, j])
    print("# conduction divergence is INDEX-SPACE, -(f1[i+1] - f1[i]), no 1/dx and no"
          " area factor")
    print("%4s %12s %12s %14s %14s %14s"
          % ("i", "T[K]", "rho", "rt_src", "rt_de", "-d(cond_f1)"))
    f1 = d["cond_f1"][k, j]
    for i in range(is_, ie + 1):
        div = -(f1[i + 1] - f1[i]) if i + 1 < f1.size else float("nan")
        print("%4d %12.5e %12.5e %14.6e %14.6e %14.6e"
              % (i, d["rt_T"][k, j, i], d["w_dens"][k, j, i],
                 d["rt_src"][k, j, i], d["rt_de"][k, j, i], div))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
