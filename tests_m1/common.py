"""Shared helpers for the <rad_m1> verification gates (T1-T6).

Everything here works in CODE UNITS.  The speed of light ``c`` and the
radiation constant ``a_r`` are never guessed from a dump: they are passed on
the command line and must match the ``<rad_m1>`` / ``<units>`` block of the
athinput that produced the dump (design note section 1).

Dump loading goes through ``vis/python/bin_convert.py`` (``read_binary`` /
``read_all_ranks_binary``).  For the ``--selftest`` modes the same loader also
accepts a ``.npz`` file with keys ``time``, ``x1``, ``x2``, ``x3`` and one
entry per variable holding a ``(nx3, nx2, nx1)`` array; that is how the
synthetic analytic data are fed to the very same analysis path.
"""

import argparse
import os
import re
import sys

import numpy as np

C_CGS = 2.99792458e10
AR_CGS = 7.5657332502e-15

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULT_BINDIR = os.path.join(_REPO, "vis", "python")


def _import_bin_convert(bin_convert_dir=None):
    """Import vis/python/bin_convert.py and return the module."""
    path = bin_convert_dir or _DEFAULT_BINDIR
    if path not in sys.path:
        sys.path.insert(0, path)
    import bin_convert  # noqa: E402
    return bin_convert


class Dump(object):
    """A uniform-mesh snapshot: cell centres plus (nx3, nx2, nx1) arrays."""

    def __init__(self, time, x1, x2, x3, data, source=""):
        self.time = float(time)
        self.x1 = np.asarray(x1, dtype=float)
        self.x2 = np.asarray(x2, dtype=float)
        self.x3 = np.asarray(x3, dtype=float)
        self.data = data
        self.source = source

    @property
    def dx1(self):
        return float(self.x1[1] - self.x1[0]) if self.x1.size > 1 else 1.0

    @property
    def dx2(self):
        return float(self.x2[1] - self.x2[0]) if self.x2.size > 1 else 1.0

    @property
    def dx3(self):
        return float(self.x3[1] - self.x3[0]) if self.x3.size > 1 else 1.0

    def has(self, name):
        return name in self.data

    def var(self, name):
        if name not in self.data:
            raise KeyError("variable '%s' not in %s (have: %s)"
                           % (name, self.source, ", ".join(sorted(self.data))))
        return self.data[name]


def centres(lo, hi, n):
    """Cell centres of a uniform grid of n cells on [lo, hi]."""
    dx = (hi - lo) / n
    return lo + dx * (np.arange(n) + 0.5)


_centres = centres


def _assemble(fd):
    """Stitch the MeshBlocks of a single-level bin dump into global arrays."""
    lev = np.asarray(fd["mb_logical"])[:, 3]
    if lev.size and (lev != lev[0]).any():
        raise ValueError("dump is multi-level (SMR/AMR); these gates are "
                         "uniform-mesh only")
    nx1, nx2, nx3 = fd["Nx1"], fd["Nx2"], fd["Nx3"]
    out = {}
    for name in fd["var_names"]:
        arr = np.zeros((nx3, nx2, nx1), dtype=float)
        for mb in range(fd["n_mbs"]):
            il, iu, jl, ju, kl, ku = fd["mb_index"][mb]
            arr[kl:ku + 1, jl:ju + 1, il:iu + 1] = fd["mb_data"][name][mb]
        out[name] = arr
    x1 = _centres(fd["x1min"], fd["x1max"], nx1)
    x2 = _centres(fd["x2min"], fd["x2max"], nx2)
    x3 = _centres(fd["x3min"], fd["x3max"], nx3)
    return x1, x2, x3, out


def load_dump(path, bin_convert_dir=None, all_ranks=False):
    """Load an AthenaK ``file_type = bin`` dump (or a synthetic .npz)."""
    if path.endswith(".npz"):
        z = np.load(path)
        skip = ("time", "x1", "x2", "x3")
        data = {k: np.asarray(z[k], dtype=float) for k in z.files
                if k not in skip}
        return Dump(float(z["time"]), z["x1"], z["x2"], z["x3"], data, path)
    bc = _import_bin_convert(bin_convert_dir)
    if all_ranks:
        fd = bc.read_all_ranks_binary(path)
    else:
        fd = bc.read_binary(path)
    x1, x2, x3, data = _assemble(fd)
    return Dump(fd["time"], x1, x2, x3, data, path)


def load_series(paths, bin_convert_dir=None, all_ranks=False):
    """Load several dumps and sort them by time."""
    dumps = [load_dump(p, bin_convert_dir, all_ranks) for p in paths]
    dumps.sort(key=lambda d: d.time)
    return dumps


def write_npz_dump(path, time, x1, x2, x3, **variables):
    """Write a synthetic dump readable by :func:`load_dump` (selftests)."""
    np.savez(path, time=np.array(float(time)), x1=np.asarray(x1),
             x2=np.asarray(x2), x3=np.asarray(x3),
             **{k: np.asarray(v, dtype=float) for k, v in variables.items()})
    return path if path.endswith(".npz") else path + ".npz"


def extract_1d(dump, name, axis=1, reduce="mid"):
    """Return (x, profile) along x1|x2|x3 from a 1-, 2- or 3-D dump.

    ``reduce`` is 'mid' (slice through the middle of the other directions)
    or 'mean' (average over them).
    """
    arr = dump.var(name)
    coord = {1: dump.x1, 2: dump.x2, 3: dump.x3}[axis]
    kax = {1: 2, 2: 1, 3: 0}[axis]
    other = [a for a in (0, 1, 2) if a != kax]
    if reduce == "mean":
        prof = arr.mean(axis=tuple(other))
    elif reduce == "mid":
        prof = arr
        for a in sorted(other, reverse=True):
            prof = np.take(prof, prof.shape[a] // 2, axis=a)
    else:
        raise ValueError("reduce must be 'mid' or 'mean'")
    return coord, np.asarray(prof, dtype=float)


def sample_2d(x1, x2, arr, px, py):
    """Bilinear sample of a (nx2, nx1) array at points (px, py)."""
    px = np.atleast_1d(np.asarray(px, dtype=float))
    py = np.atleast_1d(np.asarray(py, dtype=float))
    dx1 = x1[1] - x1[0]
    dx2 = x2[1] - x2[0]
    fi = np.clip((px - x1[0]) / dx1, 0.0, x1.size - 1.0000001)
    fj = np.clip((py - x2[0]) / dx2, 0.0, x2.size - 1.0000001)
    i0 = fi.astype(int)
    j0 = fj.astype(int)
    ti = fi - i0
    tj = fj - j0
    i1 = np.minimum(i0 + 1, x1.size - 1)
    j1 = np.minimum(j0 + 1, x2.size - 1)
    return ((1 - ti) * (1 - tj) * arr[j0, i0] + ti * (1 - tj) * arr[j0, i1]
            + (1 - ti) * tj * arr[j1, i0] + ti * tj * arr[j1, i1])


def l1_rel(a, b):
    """Relative L1 norm ||a - b||_1 / ||b||_1 (unweighted, equal cells)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    den = np.abs(b).sum()
    if den == 0.0:
        return np.abs(a - b).sum()
    return np.abs(a - b).sum() / den


def l1_abs(a, b, dx=1.0):
    return float(np.abs(np.asarray(a) - np.asarray(b)).sum() * dx)


def moments(x, w, background=0.0):
    """Weighted mean and variance of ``x`` with weights ``w - background``."""
    wp = np.asarray(w, dtype=float) - background
    wp = np.maximum(wp, 0.0)
    tot = wp.sum()
    if tot <= 0.0:
        raise ValueError("no positive weight left after background removal")
    mean = float((wp * x).sum() / tot)
    var = float((wp * (x - mean) ** 2).sum() / tot)
    return mean, var


def linfit(t, y):
    """Least-squares slope and intercept of y(t)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    a = np.vstack([t, np.ones_like(t)]).T
    sol, _, _, _ = np.linalg.lstsq(a, y, rcond=None)
    return float(sol[0]), float(sol[1])


def nyquist_amplitude(prof):
    """Signed projection of a profile on the odd-even mode (-1)**i."""
    prof = np.asarray(prof, dtype=float)
    sgn = np.where(np.arange(prof.size) % 2 == 0, 1.0, -1.0)
    return float((prof * sgn).mean())


def nyquist_power_fraction(prof):
    """Power in the Nyquist mode / total fluctuation power."""
    prof = np.asarray(prof, dtype=float)
    f = np.fft.rfft(prof - prof.mean())
    p = np.abs(f) ** 2
    tot = p.sum()
    if tot <= 0.0:
        return 0.0
    return float(p[-1] / tot)


def read_history(path):
    """Read an AthenaK history file.  Returns (names, array[nrow, ncol])."""
    names = None
    rows = []
    with open(path) as fp:
        for line in fp:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                found = re.findall(r"\[(\d+)\]=(\S+)", s)
                if found:
                    names = [n for _, n in sorted(found,
                                                  key=lambda kv: int(kv[0]))]
                continue
            rows.append([float(v) for v in s.split()])
    arr = np.array(rows, dtype=float)
    if names is None:
        names = ["col%d" % (i + 1) for i in range(arr.shape[1])]
    return names, arr


def hist_column(names, arr, key):
    """Column of a history file by name (substring match) or 1-based index."""
    if isinstance(key, int) or (isinstance(key, str) and key.isdigit()):
        return arr[:, int(key) - 1]
    if key in names:
        return arr[:, names.index(key)]
    hits = [i for i, n in enumerate(names) if key in n]
    if len(hits) == 1:
        return arr[:, hits[0]]
    raise KeyError("history column '%s' not found (have: %s)"
                   % (key, ", ".join(names)))


def base_parser(description):
    """Argparse boilerplate shared by every gate script."""
    p = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--bin-convert-dir", default=None,
                   help="directory holding bin_convert.py "
                        "(default: <repo>/vis/python)")
    p.add_argument("--all-ranks", action="store_true",
                   help="dumps are one-file-per-rank; pass the rank0 file")
    p.add_argument("--selftest", action="store_true",
                   help="run on synthetic analytic data instead of dumps")
    p.add_argument("--quiet", action="store_true",
                   help="print only the verdict line")
    return p


def report(args, text):
    if not getattr(args, "quiet", False):
        print(text)


def verdict(ok, text):
    """Print the one-line verdict and exit (0 on PASS, 1 on FAIL)."""
    print(("PASS " if ok else "FAIL ") + text)
    sys.exit(0 if ok else 1)
