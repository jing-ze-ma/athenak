"""rsgE page: plotting helpers on the NATIVE cubed-sphere cells (no resampling of the data).
ring_xy   : pcolormesh of a native ring (angle edges x radial faces) in the plane, every cell a true polar quad.
ring_unwrap: same ring unwrapped, angle [deg] on x and r on y (photosphere band).
moll      : Mollweide image of a (6, N, N) shell where every screen pixel is inverse-projected and takes the value of the
            cubed-sphere cell it falls in (so each cell is drawn with its true outline at screen resolution)."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LogNorm, Normalize
from matplotlib.patches import Ellipse

FR = {0: ((0, 1, 0), (0, 0, 1), (1, 0, 0)), 1: ((-1, 0, 0), (0, 0, 1), (0, 1, 0)), 2: ((0, -1, 0), (0, 0, 1), (-1, 0, 0)),
      3: ((0, 1, 0), (-1, 0, 0), (0, 0, 1)), 4: ((1, 0, 0), (0, 0, 1), (0, -1, 0)), 5: ((0, 1, 0), (1, 0, 0), (0, 0, -1))}
FA = np.array([FR[p][0] for p in range(6)], float)
FB = np.array([FR[p][1] for p in range(6)], float)
FN = np.array([FR[p][2] for p in range(6)], float)
FLOORC = '#8a8f96'          # T_gas floor cells


def cmap_floor(name):
    c = plt.get_cmap(name).copy()
    c.set_under(FLOORC)
    return c


def norm_of(kind, lim):
    if kind == 'sym':
        return TwoSlopeNorm(0, -lim, lim)
    if kind == 'log':
        return LogNorm(lim[0], lim[1])
    return Normalize(lim[0], lim[1])


def under_value(norm):
    if isinstance(norm, LogNorm):
        return norm.vmin * 1e-3
    return norm.vmin - 10 * (abs(norm.vmax - norm.vmin) + 1)


def with_floor(q, fl, norm):
    """clip into the colour range (so only floor cells use the 'under' colour), then mark floor cells"""
    if fl is None:
        return q
    lo, hi = norm.vmin, norm.vmax
    q = np.clip(q, lo * (1 + 1e-6) if isinstance(norm, LogNorm) else lo + 1e-9 * abs(hi - lo), hi)
    return np.where(fl, under_value(norm), q)


def shell_mean(sh):
    """solid-angle weighted mean of a (6, N, N) native shell"""
    N = sh.shape[1]
    x = np.tan((-1 + (np.arange(N) + 0.5) * 2.0 / N) * np.pi / 4)
    X, Y = np.meshgrid(x, x, indexing='xy')
    w = (1 + X ** 2) * (1 + Y ** 2) / (1 + X ** 2 + Y ** 2) ** 1.5
    return float((np.asarray(sh, float) * w[None]).sum() / (6 * w.sum()))


def ring_mesh(ang, rf, i0, i1):
    AA, RR = np.meshgrid(ang, rf[i0:i1 + 1], indexing='ij')
    return RR * np.cos(AA), RR * np.sin(AA)


def ring_xy(ax, ang, rf, q, i0, i1, cmap, norm, swap=False):
    X, Y = ring_mesh(ang, rf, i0, i1)
    return ax.pcolormesh(X, Y, q[:, i0:i1], cmap=cmap, norm=norm, shading='flat', rasterized=True, linewidth=0,
                         antialiased=False)


def ring_unwrap(ax, ang, rf, q, i0, i1, cmap, norm):
    a = np.degrees(ang)
    return ax.pcolormesh(a, rf[i0:i1 + 1], q[:, i0:i1].T, cmap=cmap, norm=norm, shading='flat', rasterized=True,
                         linewidth=0, antialiased=False)


def rrange(rf, lo, hi):
    """cell index range [i0, i1) with cells overlapping lo..hi (rf in R_ph units)"""
    i0 = max(0, int(np.searchsorted(rf, lo, side='right')) - 1)
    i1 = min(len(rf) - 1, int(np.searchsorted(rf, hi, side='left')))
    return i0, i1


_LUT = {}


def moll_lut(N, W, H, cache=None):
    """pixel -> flat cell index (p*N*N + k*N + j) for a W x H Mollweide image; -1 outside the ellipse"""
    key = (N, W, H)
    if key in _LUT:
        return _LUT[key]
    fn = os.path.join(cache, 'moll_lut_%d_%dx%d.npy' % key) if cache else None
    if fn and os.path.exists(fn):
        _LUT[key] = np.load(fn)
        return _LUT[key]
    s2 = np.sqrt(2.)
    x = (np.arange(W) + 0.5) / W * 4 * s2 - 2 * s2
    y = s2 - (np.arange(H) + 0.5) / H * 2 * s2
    X, Y = np.meshgrid(x, y)
    inside = (X / (2 * s2)) ** 2 + (Y / s2) ** 2 < 1
    th = np.arcsin(np.clip(Y / s2, -1, 1))
    lat = np.arcsin(np.clip((2 * th + np.sin(2 * th)) / np.pi, -1, 1))
    lon = np.pi * X / (2 * s2 * np.maximum(np.cos(th), 1e-12))
    inside &= np.abs(lon) <= np.pi
    q = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)
    p = np.argmax(q @ FN.T, axis=-1)
    qa = np.einsum('...c,...c->...', q, FA[p]); qb = np.einsum('...c,...c->...', q, FB[p]); qn = np.einsum('...c,...c->...', q, FN[p])
    j = np.clip(((np.arctan(qa / qn) / (np.pi / 4) + 1) * 0.5 * N).astype(int), 0, N - 1)
    k = np.clip(((np.arctan(qb / qn) / (np.pi / 4) + 1) * 0.5 * N).astype(int), 0, N - 1)
    lut = np.where(inside, p * N * N + k * N + j, -1).astype(np.int32)
    if fn:
        np.save(fn, lut)
    _LUT[key] = lut
    return lut


def moll(ax, sh, cmap, norm, W=1000, cache=None, floor=None):
    """draw (6,N,N) shell sh on a plain axes as a native-cell Mollweide image; returns a ScalarMappable"""
    N = sh.shape[1]
    H = W // 2
    lut = moll_lut(N, W, H, cache)
    flat = np.asarray(sh, float).ravel()
    img = np.where(lut >= 0, flat[np.maximum(lut, 0)], np.nan)
    if floor is not None:
        fl = np.asarray(floor).ravel()[np.maximum(lut, 0)] & (lut >= 0)
        img = with_floor(img, fl, norm)                # floor cells -> cmap 'under' colour
    s2 = np.sqrt(2.)
    m = ax.imshow(np.ma.masked_invalid(img), cmap=cmap, norm=norm, extent=(-2 * s2, 2 * s2, -s2, s2),
                  interpolation='nearest', origin='upper')
    ax.add_patch(Ellipse((0, 0), 4 * s2, 2 * s2, fill=False, ec='#6b7480', lw=0.6))
    for lo in np.radians([-135, -90, -45, 0, 45, 90, 135]):              # graticule every 45 deg (cube edges at +-45, +-135)
        la = np.linspace(-np.pi / 2, np.pi / 2, 181)
        th = _theta(la)
        ax.plot(2 * s2 / np.pi * lo * np.cos(th), s2 * np.sin(th), color='#6b7480', lw=0.35)
    for la0 in np.radians([-60, -30, 0, 30, 60]):
        th = _theta(np.array([la0]))[0]
        ax.plot([-2 * s2 * np.cos(th), 2 * s2 * np.cos(th)], [s2 * np.sin(th)] * 2, color='#6b7480', lw=0.35)
    ax.set_xlim(-2 * s2 * 1.01, 2 * s2 * 1.01); ax.set_ylim(-s2 * 1.02, s2 * 1.02)
    ax.set_aspect('equal'); ax.axis('off')
    return m


def _theta(lat):
    th = lat.copy()
    for _ in range(30):
        f = 2 * th + np.sin(2 * th) - np.pi * np.sin(lat)
        th = th - f / np.maximum(2 + 2 * np.cos(2 * th), 1e-9)
    return th
