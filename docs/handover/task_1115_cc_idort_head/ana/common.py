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


def _merge_axis(values, tol=1e-9):
    """Sorted unique coordinate values, merging duplicates within a tolerance.

    ``values`` is the concatenation of the cell centres of every MeshBlock
    along one axis, so the same centre appears once per block that owns it.
    Blocks of one level and one size produce bit-identical centres, but the
    merge is done with a tolerance anyway so that a dump written in single
    precision (which every ``file_type = bin`` dump is) still collapses.
    """
    v = np.sort(np.asarray(values, dtype=float))
    if v.size == 0:
        return v
    span = max(float(v[-1] - v[0]), abs(float(v[0])), 1.0)
    keep = [v[0]]
    for a in v[1:]:
        if a - keep[-1] > tol * span:
            keep.append(a)
    return np.array(keep, dtype=float)


def _assemble(fd):
    """Stitch the MeshBlocks of a single-level bin dump into global arrays.

    ``mb_index`` is NOT a global index range: it is the index range of the
    block's output slab inside the block itself (it is ``ois..oie`` etc., and
    for a collapsed or sliced direction it comes back as ``-2``, which is why
    the naive version of this function wrote nothing at all for a 1-D or 2-D
    dump and every gate saw zeros).  The block position is therefore taken
    from ``mb_geometry`` (the block's own x1min/x1max/... ), which is correct
    for 1-, 2- and 3-D, for many blocks, and for sliced outputs.
    """
    lev = np.asarray(fd["mb_logical"])[:, 3]
    if lev.size and (lev != lev[0]).any():
        raise ValueError("dump is multi-level (SMR/AMR); these gates are "
                         "uniform-mesh only")
    names = list(fd["var_names"])
    nmb = int(fd["n_mbs"])
    geom = np.asarray(fd["mb_geometry"], dtype=float)
    shapes = [np.asarray(fd["mb_data"][names[0]][mb]).shape
              for mb in range(nmb)]
    # per-block cell centres along each axis
    cen = [[], [], []]
    for mb in range(nmb):
        n3b, n2b, n1b = shapes[mb]
        cen[0].append(_centres(geom[mb, 0], geom[mb, 1], n1b))
        cen[1].append(_centres(geom[mb, 2], geom[mb, 3], n2b))
        cen[2].append(_centres(geom[mb, 4], geom[mb, 5], n3b))
    x1 = _merge_axis(np.concatenate(cen[0]))
    x2 = _merge_axis(np.concatenate(cen[1]))
    x3 = _merge_axis(np.concatenate(cen[2]))
    nx1, nx2, nx3 = x1.size, x2.size, x3.size

    offs = []
    for mb in range(nmb):
        o1 = int(np.abs(x1 - cen[0][mb][0]).argmin())
        o2 = int(np.abs(x2 - cen[1][mb][0]).argmin())
        o3 = int(np.abs(x3 - cen[2][mb][0]).argmin())
        n3b, n2b, n1b = shapes[mb]
        if o1 + n1b > nx1 or o2 + n2b > nx2 or o3 + n3b > nx3:
            raise ValueError("MeshBlock %d does not fit the assembled grid "
                             "(overlapping or non-uniform blocks?)" % mb)
        offs.append((o3, o2, o1))

    out = {}
    for name in names:
        arr = np.full((nx3, nx2, nx1), np.nan, dtype=float)
        for mb in range(nmb):
            o3, o2, o1 = offs[mb]
            n3b, n2b, n1b = shapes[mb]
            arr[o3:o3 + n3b, o2:o2 + n2b, o1:o1 + n1b] = \
                np.asarray(fd["mb_data"][name][mb], dtype=float)
        if np.isnan(arr).any():
            raise ValueError("the MeshBlocks of this dump do not tile the "
                             "domain (%d of %d cells unfilled)"
                             % (int(np.isnan(arr).sum()), arr.size))
        out[name] = arr
    return x1, x2, x3, out


_TAB_COORDS = ("x1v", "x2v", "x3v")


def load_tab(path):
    """Load an AthenaK ``file_type = tab`` 1-D (or sliced) output.

    The bin writer is SINGLE precision, so any gate that needs more than
    ~7 digits (the T4b gas-temperature drift, the T4 advected-vs-static
    difference) runs with ``file_type = tab`` and a wide ``data_format``
    instead.  Two traps in that writer (``formatted_table.cpp``): the header
    line names a coordinate pair for every direction that is not *sliced*,
    while the data rows carry one only for the directions that actually have
    more than one cell -- so for a 1-D run the header has six coordinate
    tokens and the rows two; and with several MeshBlocks or ranks the rows
    come out in block order, not in coordinate order.  Both are handled here.
    """
    time = 0.0
    header = None
    rows = []
    with open(path) as fp:
        for line in fp:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                mt = re.search(r"time=(\S+)", s)
                if mt:
                    time = float(mt.group(1))
                elif header is None:
                    header = s.lstrip("#").split()
                continue
            rows.append([float(v) for v in s.split()])
    if header is None or not rows:
        raise ValueError("%s is not an AthenaK tab file" % path)
    dat = np.array(rows, dtype=float)
    ncol = dat.shape[1]
    # header: gid [i x1v] [j x2v] [k x3v] var...
    hcoord = [n for n in header if n in _TAB_COORDS]
    labels = header[header.index(hcoord[-1]) + 1:] if hcoord else header[1:]
    ndim = (ncol - len(labels) - 1) // 2
    if ndim < 1 or 1 + 2 * ndim + len(labels) != ncol:
        raise ValueError("%s: cannot match %d columns to the header %s"
                         % (path, ncol, " ".join(header)))
    axes = hcoord[:ndim]
    coords = {}
    idx = {}
    for d in range(ndim):
        col = dat[:, 2 + 2 * d]
        coords[axes[d]] = _merge_axis(col)
        idx[axes[d]] = np.abs(coords[axes[d]][None, :]
                              - col[:, None]).argmin(axis=1)
    x1 = coords.get("x1v", np.array([0.0]))
    x2 = coords.get("x2v", np.array([0.0]))
    x3 = coords.get("x3v", np.array([0.0]))
    i1 = idx.get("x1v", np.zeros(dat.shape[0], dtype=int))
    i2 = idx.get("x2v", np.zeros(dat.shape[0], dtype=int))
    i3 = idx.get("x3v", np.zeros(dat.shape[0], dtype=int))
    data = {}
    for n, lab in enumerate(labels):
        arr = np.full((x3.size, x2.size, x1.size), np.nan)
        arr[i3, i2, i1] = dat[:, 1 + 2 * ndim + n]
        if np.isnan(arr).any():
            raise ValueError("%s: %d cells missing from the table"
                             % (path, int(np.isnan(arr).sum())))
        data[lab] = arr
    return Dump(time, x1, x2, x3, data, path)


def load_dump(path, bin_convert_dir=None, all_ranks=False):
    """Load an AthenaK ``file_type = bin`` dump (or a synthetic .npz)."""
    if path.endswith(".tab"):
        return load_tab(path)
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


# ---------------------------------------------------------------------------
# loader selftest: reads REAL dumps produced by the code (data/), not
# synthetic arrays, because what was broken was the reading of a real file.
# data/loader_1d.{bin,tab} and data/loader_2d.bin are the t = 0 output of
# data/loader_1d.athinput (problem/m1_test = thick_pulse, a Gaussian in E on
# top of a background), 32 cells in 2 MeshBlocks and 16 x 8 cells in 4.

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _loader_reference(x1):
    """The initial condition of data/loader_1d.athinput, analytically."""
    return 1.0e-3 + 1.0 * np.exp(-((x1 - 0.5) / 0.05) ** 2)


def loader_selftest(quiet=False, fail=False):
    """Check the dump loader against the analytic IC.  Returns True on PASS."""
    ok = True
    msgs = []

    d1 = load_dump(os.path.join(_DATA, "loader_1d.bin"))
    ref = _loader_reference(d1.x1)
    e = d1.var("m1_e")
    err = float(np.abs(e[0, 0, :] - ref).max() / ref.max())
    if fail:
        err *= 1.0e6
    msgs.append("1-D bin: shape=%s nx1=%d max|E-E_exact|/max(E)=%.2e"
                % (e.shape, d1.x1.size, err))
    ok = ok and e.shape == (1, 1, 32) and err < 1e-6

    dt = load_dump(os.path.join(_DATA, "loader_1d.tab"))
    et = dt.var("m1_e")
    dbt = float(np.abs(et - e).max())
    # the tab file is written at %24.16e, the bin file in single precision
    errt = float(np.abs(et[0, 0, :] - ref).max() / ref.max())
    if fail:
        errt *= 1.0e6
    msgs.append("1-D tab: nx1=%d max|E_tab-E_exact|/max(E)=%.2e  "
                "max|tab-bin|=%.2e (single-precision bin)"
                % (dt.x1.size, errt, dbt))
    ok = ok and et.shape == (1, 1, 32) and errt < 1e-14 and dbt < 1e-6

    d2 = load_dump(os.path.join(_DATA, "loader_2d.bin"))
    e2 = d2.var("m1_e")
    ref2 = _loader_reference(d2.x1)
    err2 = float(np.abs(e2[0] - ref2[None, :]).max() / ref2.max())
    if fail:
        err2 *= 1.0e6
    msgs.append("2-D bin: shape=%s (4 MeshBlocks)  "
                "max|E-E_exact|/max(E)=%.2e" % (e2.shape, err2))
    ok = ok and e2.shape == (1, 8, 16) and err2 < 1e-6

    # nothing may come back identically zero: that was the symptom of the bug
    for nm, dd in (("1-D bin", d1), ("1-D tab", dt), ("2-D bin", d2)):
        if not np.any(dd.var("m1_e") != 0.0):
            ok = False
            msgs.append("%s: every cell is zero" % nm)

    if not quiet:
        for m in msgs:
            print(m)
    return ok


if __name__ == "__main__":
    _ap = argparse.ArgumentParser(description="common.py loader selftest")
    _ap.add_argument("--selftest", action="store_true")
    _ap.add_argument("--selftest-fail", action="store_true")
    _ap.add_argument("--quiet", action="store_true")
    _a = _ap.parse_args()
    _ok = loader_selftest(_a.quiet, _a.selftest_fail)
    verdict(_ok, "common.load_dump: 1-D/2-D, multi-block, bin and tab")
