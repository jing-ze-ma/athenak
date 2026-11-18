#!/usr/bin/env python3
"""accsp_prod_page (11-02): per-dump reduction of the sp full-sphere accretor PRODUCTION (accsp_1102/prod/fresh; READ ONLY on the run).
usage: reduce_sp.py <out_dir> <hydro_w.bin> <m1.bin|-> <tag>   -> out_dir/red_<tag>.npz   (tag = a_NNNNN static, b_NNNNN stream)
Same quantities as acc3d_cs128_page/src/reduce_cs.py, on the native sp grid (8 MeshBlocks along phi, 64 x 116 x 500 each):
velx = v_r, vely = v_theta, velz = v_phi (rotating frame). Orbital plane = the two theta rows next to the equator averaged;
x-z plane = phi 0 (cells k 511, 0) and phi pi (k 255, 256) half planes; log T_gas from the run's general EOS (eos_gen numpy port,
as the cs page); log P_rad/P_gas (E/3) when an m1 dump has the same time; Rosseland photosphere R(tau=2/3) per (phi, theta)
column (MESA kap table of the deck, cells r > RTRC); L1 tip column (phi 0, equator, 2 x 2 cells); floor-class counts
(rho < 1e-7 code); KE by |latitude| band x r bin from accsp_1102/tools/keband.py (called as is); phi-rms of v_phi in the
polar theta rows at r ~ 6.8 (the ring-averaged rows, the 11-02 lateral-KE source)."""
import os
import re
import subprocess
import sys
import time
import types
import numpy as np
J = '/resnick/groups/carnegie_poc/jingze'
sys.modules.setdefault('h5py', types.ModuleType('h5py'))
sys.path.insert(0, J + '/rsgE_page/src')
sys.path.insert(0, J + '/wt_accsp/vis/python')
import bin_convert as bc  # noqa: E402
import eos_gen  # noqa: E402

VER = 2   # trc-1102: + map_s0ph (s0 at the photosphere cell), prof_m / prof_ms0 / prof_ms1 (shell masses), eq/me s1
OUT, FH, FM, TAG = sys.argv[1], sys.argv[2], (sys.argv[3] if sys.argv[3] != '-' else None), sys.argv[4]
RHOU = 5.1099127e-08
EU = RHOU * 1e10
LU = 6.957e10
RTRC = 8.4
KEBAND = J + '/accsp_1102/tools/keband.py'
OPF = J + '/acc3d_run12_move/in/rosseland_mesakap_gs98_x0.70_z0.02.txt'
hdr = [ln for ln in open(OPF) if ln.startswith('#')]
_g = [ln for ln in hdr if ln.startswith('# grid')][0]
nT, nD, lT0, dlT, lD0, dlD = [float(x) for x in hdr[hdr.index(_g) + 1][1:].split()]
nT, nD = int(nT), int(nD)
KAP = np.loadtxt(OPF, comments='#').reshape(nT, nD)


def kR(rho, T):
    i = (np.log10(T) - lT0) / dlT; j = (np.log10(rho) - lD0) / dlD
    i0 = np.clip(np.floor(i).astype(int), 0, nT - 2); j0 = np.clip(np.floor(j).astype(int), 0, nD - 2)
    fi = np.clip(i - i0, 0, 1); fj = np.clip(j - j0, 0, 1)
    return 10 ** ((1 - fi) * (1 - fj) * KAP[i0, j0] + fi * (1 - fj) * KAP[i0 + 1, j0] + (1 - fi) * fj * KAP[i0, j0 + 1] + fi * fj * KAP[i0 + 1, j0 + 1])


def assemble(d, var):
    """-> (nphi, ntheta, nr), MeshBlocks ordered by their x3 position"""
    order = np.argsort([float(d['mb_x3f'][m][0]) for m in range(len(d['mb_logical']))])
    return np.concatenate([np.asarray(d['mb_data'][var][m], np.float32) for m in order], axis=0)


def keband(fn):
    """run accsp_1102/tools/keband.py sp and parse -> (2, 5) [radial row sums | lateral by band summed over r] and lateral (5, 4)"""
    o = subprocess.run([sys.executable, KEBAND, 'sp', fn], capture_output=True, text=True).stdout
    lat = []; rad = []; tot = None
    for ln in o.splitlines():
        m = re.match(r'\s*(\d+)-\s*(\d+) deg\s+(.*)\|\s*([0-9.eE+-]+)', ln)
        if m:
            lat.append([float(x) for x in m.group(3).split()]); rad.append(float(m.group(4)))
        m = re.search(r'KE total ([0-9.eE+-]+) radial ([0-9.eE+-]+) lateral ([0-9.eE+-]+)', ln)
        if m:
            tot = [float(m.group(i)) for i in (1, 2, 3)]
    return np.array(lat), np.array(rad), np.array(tot if tot else [np.nan] * 3)


t0 = time.time()
h = bc.read_binary(FH)
t = float(h['time'])
r = np.asarray(h['mb_x1v'][0], float); rf = np.asarray(h['mb_x1f'][0], float)
thf = np.asarray(h['mb_x2f'][0], float); th = np.asarray(h['mb_x2v'][0], float)
phf = np.sort(np.unique(np.concatenate([np.asarray(h['mb_x3f'][m], float) for m in range(len(h['mb_logical']))]).round(12)))
rho = assemble(h, 'dens'); vr = assemble(h, 'velx'); vt = assemble(h, 'vely'); vp = assemble(h, 'velz'); ei = assemble(h, 'eint')
s0 = assemble(h, 's_00') if 's_00' in h['var_names'] else None
s1 = assemble(h, 's_01') if 's_01' in h['var_names'] else None
del h
NP, NT, NR = rho.shape
E = None
if FM:
    m = bc.read_binary(FM)
    if abs(float(m['time']) - t) > 1e-6:
        print('m1 time mismatch', FM, m['time'], t); FM = None
    else:
        E = assemble(m, 'm1_e')
    del m
TBL = eos_gen.Table(logd=(-18.0, 0.0), logt=(3.5, 7.6), dlog=0.02, xh=0.70, yhe=0.28)
je = (NT // 2 - 1, NT // 2)
k0, kp = (NP - 1, 0), (NP // 2 - 1, NP // 2)


def ring_eq(q):
    return 0.5 * (q[:, je[0], :] + q[:, je[1], :])                     # (nphi, nr)


def ring_me(q):
    return np.stack([0.5 * (q[k0[0]] + q[k0[1]]), 0.5 * (q[kp[0]] + q[kp[1]])])   # (2 [phi 0, phi pi], ntheta, nr)


def tgas(rr, ee):
    T, _ = TBL.temperature(rr.astype(float) * RHOU, ee.astype(float) * EU)
    return np.asarray(T, float)


out = dict(ver=VER, t=t, r=r, rf=rf, thf=thf, th=th, phf=phf, has_m1=E is not None, idx=TAG, src=FH)
for key, ring in (('eq', ring_eq), ('me', ring_me)):
    rr = ring(rho); ee = ring(ei); T = tgas(rr, ee)
    out[key + '_lrho'] = np.log10(np.maximum(rr, 1e-30)).astype(np.float32)
    out[key + '_lT'] = np.log10(np.maximum(T, 1.0)).astype(np.float32)
    out[key + '_vr'] = ring(vr).astype(np.float32)
    out[key + '_vphi'] = ring(vp).astype(np.float32)
    if s0 is not None:
        out[key + '_s0'] = ring(s0).astype(np.float32)
    if s1 is not None:
        out[key + '_s1'] = ring(s1).astype(np.float32)
    if E is not None:
        pg = np.asarray(TBL.pressure(rr.astype(float) * RHOU, T), float) * rr * RHOU
        out[key + '_lbeta'] = np.log10(np.maximum(ring(E).astype(float) * EU / 3.0, 1e-30) / np.maximum(pg, 1e-30)).astype(np.float32)
# photosphere per (phi, theta) column, tau inward from x1max, cells r > RTRC
i0 = int(np.searchsorted(r, RTRC))
Rph = np.empty((NP, NT), np.float32); Tph = np.empty((NP, NT), np.float32)
S0ph = np.full((NP, NT), np.nan, np.float32)
drc = np.diff(rf)[i0:] * LU; rfl = rf[i0:]
for k in range(0, NP, 64):
    rr = rho[k:k + 64, :, i0:].reshape(-1, NR - i0).astype(float) * RHOU
    ee = ei[k:k + 64, :, i0:].reshape(-1, NR - i0).astype(float) * EU
    T, _ = TBL.temperature(rr, ee); T = np.asarray(T, float).reshape(rr.shape)
    dtau = kR(rr, np.maximum(T, 10 ** lT0)) * rr * drc[None, :]
    tin = np.cumsum(dtau[:, ::-1], axis=1)[:, ::-1]
    ok = tin >= 2.0 / 3.0
    n = (NR - i0) - 1 - np.argmax(ok[:, ::-1], axis=1)
    c = np.arange(len(n)); tout = tin[c, n] - dtau[c, n]
    R = rfl[n + 1] - (2.0 / 3.0 - tout) / np.maximum(dtau[c, n], 1e-300) * (rfl[n + 1] - rfl[n])
    R = np.where(ok.any(axis=1), R, rfl[0])
    Rph[k:k + 64] = R.reshape(-1, NT); Tph[k:k + 64] = T[c, n].reshape(-1, NT)
    if s0 is not None:
        S0ph[k:k + 64] = s0[k:k + 64, :, i0:].reshape(-1, NR - i0)[c, n].reshape(-1, NT)
    del rr, ee, T, dtau, tin, ok
out['map_Rph'] = Rph; out['map_Tph'] = Tph
out['map_s0ph'] = S0ph
# trc-1102: shell masses [code] (rho s dV summed over phi, theta): where the stream gas sits in r
_w = np.outer(np.diff(phf), -np.diff(np.cos(thf)))[:, :, None] * (np.diff(rf ** 3) / 3.0)[None, None, :]
out['prof_m'] = np.einsum('ijk,ijk->k', rho.astype(float), _w)
if s0 is not None:
    out['prof_ms0'] = np.einsum('ijk,ijk->k', rho.astype(float) * s0, _w)
if s1 is not None:
    out['prof_ms1'] = np.einsum('ijk,ijk->k', rho.astype(float) * s1, _w)
del _w
w = np.outer(np.diff(phf), -np.diff(np.cos(thf)))                     # solid-angle weights for the percentiles
def wpct(q, p):
    s = np.argsort(q, axis=None); cw = np.cumsum(w.ravel()[s]); return float(q.ravel()[s][np.searchsorted(cw, p / 100 * cw[-1])])
out['Rph_med'] = wpct(Rph, 50); out['Rph_p2'] = wpct(Rph, 2); out['Rph_p98'] = wpct(Rph, 98); out['Tph_med'] = wpct(Tph, 50)
# L1 tip: phi 0, equator, 2 x 2 cells
out['tip_rho'] = rho[[k0[0], k0[1]]][:, [je[0], je[1]], :].mean(axis=(0, 1)).astype(np.float32)
out['tip_vr'] = vr[[k0[0], k0[1]]][:, [je[0], je[1]], :].mean(axis=(0, 1)).astype(np.float32)
# floor-class gas
fl = rho < 1e-7
out['n_floor'] = int(fl.sum()); out['n_cells'] = int(rho.size)
eor = np.where(fl, ei / np.maximum(rho, 1e-30), 0)
out['n_floor_hot'] = int((eor > 1e4).sum()); out['floor_eor_max'] = float(eor.max())
band = [(0, 9.0), (9.0, 9.8), (9.8, 13.9), (13.9, 16.3), (16.3, 99)]
out['floor_band'] = np.array([int(fl[..., (r >= lo) & (r < hi)].sum()) for lo, hi in band])
out['nonfinite'] = int((~np.isfinite(rho)).sum() + (~np.isfinite(ei)).sum() + (~np.isfinite(vr)).sum())
out['vmax'] = float(np.nanmax(np.abs(vr)))
# polar rows: phi-rms of v_phi and v_theta at r ~ 6.8 (rows 0..9 from each pole, averaged over the two poles)
ip = int(np.argmin(abs(r - 6.8)))
out['pole_th'] = np.degrees(th[:10])
out['pole_vphi_rms'] = 0.5 * (vp[:, :10, ip].std(axis=0) + vp[:, ::-1, ip][:, :10].std(axis=0))
out['pole_vth_rms'] = 0.5 * (vt[:, :10, ip].std(axis=0) + vt[:, ::-1, ip][:, :10].std(axis=0))
del rho, vr, vt, vp, ei, s0, s1, E
lat, rad, tot = keband(FH)
out['keb_lat'] = lat; out['keb_rad'] = rad; out['keb_tot'] = tot
o = os.path.join(OUT, 'red_%s.npz' % TAG)
np.savez(o + '.tmp.npz', **out)
os.replace(o + '.tmp.npz', o)
print('reduced', TAG, 't %.5f' % t, 'm1' if out['has_m1'] else 'no-m1', 'Rph med %.4f p98 %.3f' % (out['Rph_med'], out['Rph_p98']),
      'floor %d' % out['n_floor'], 'KE lat %.3e' % tot[2], '%.0f s' % (time.time() - t0), flush=True)
