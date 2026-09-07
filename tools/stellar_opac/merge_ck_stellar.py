"""Splice the Rosseland mean DERIVED FROM THE CORRELATED-K TABLE into the stellar
(OPLIB+AESOPUS) one, so that the radiative-diffusion operator and the two-stream solver
use the SAME opacity where they overlap.

Why: AthenaK's tau blend hands the flux from the two-stream to diffusion between
rad_tau_lo and rad_tau_hi.  If the two operators disagree about kappa they also disagree
about tau, so the handover sits at a different physical depth for each and is not
conservative.  Measured on a red giant envelope, the two datasets differ by ~30 % above
2700 K (mostly a metallicity mismatch) and by 1-3 ORDERS OF MAGNITUDE below 1500 K, where
the stellar tables carry solid grains and the correlated-k tables carry none.

So: use the ck-derived mean where the k-table is valid and grain-free chemistry is right
(--tlo .. --thi, default 2000-5000 K), the stellar table above --thi2 (the k-table stops
at 6100 K while a stellar envelope reaches 1e6 K) and BELOW --tlo, where grains dominate
and the k-table has none, with linear blends in log T across both joins.

The ck table comes from the code itself: run the problem generator once with
problem/opac_compare=<file>, which writes log10 kappa_R on the k-table's own (log T,
log p) grid.  Converting it to the (log T, log rho) axis the conduction module reads
needs the run's own EOS, dumped with <block>/eos_table_dump.

usage: merge_ck_stellar.py CK_KR STELLAR EOS_TABLE OUT [--tlo 2000 --thi 5000 --thi2 6000]
"""
import argparse
import numpy as np


def read_ck(fn):
    rows = [ln for ln in open(fn) if not ln.startswith('#')]
    nT, nP = map(int, rows[0].split())
    v = [float(x) for x in rows[1:]]
    return np.array(v[:nT]), np.array(v[nT:nT+nP]), np.array(v[nT+nP:]).reshape(nT, nP)


def read_grid(fn):
    lines = open(fn).read().splitlines()
    h = [ln for ln in lines if ln.startswith('#')][3].split()[1:]
    nT, nD = int(h[0]), int(h[1])
    lt0, dlt, ld0, dld = (float(x) for x in h[2:])
    K = np.array([float(ln) for ln in lines
                  if not ln.startswith('#')]).reshape(nT, nD)
    return lt0 + dlt*np.arange(nT), ld0 + dld*np.arange(nD), K, (lt0, dlt, ld0, dld)


def read_eos(fn):
    lines = open(fn).read().splitlines()
    h = [ln for ln in lines if ln.startswith('#')][2].split()[1:]
    nx, ny = int(h[0]), int(h[1])
    xmin, dx, ymin, dy = (float(x) for x in h[2:])
    lpr = np.array([float(r.split()[1]) for r in lines
                    if not r.startswith('#')])
    return xmin + dx*np.arange(nx), ymin + dy*np.arange(ny), lpr.reshape(ny, nx)


ap = argparse.ArgumentParser()
ap.add_argument('ck')
ap.add_argument('stellar')
ap.add_argument('eos')
ap.add_argument('out')
ap.add_argument('--tlo', type=float, default=2000.0)
ap.add_argument('--thi', type=float, default=5000.0)
ap.add_argument('--thi2', type=float, default=6000.0)
a = ap.parse_args()

cT, cP, cK = read_ck(a.ck)
sT, sD, sK, grid = read_grid(a.stellar)
elrho, elT, elpr = read_eos(a.eos)

# log10 p on the stellar table's own (log T, log rho) nodes, from the run's EOS
LT, LD = np.meshgrid(sT, sD, indexing='ij')
jT = np.clip(((LT - elT[0])/(elT[1]-elT[0])).astype(int), 0, len(elT)-1)
iD = np.clip(((LD - elrho[0])/(elrho[1]-elrho[0])).astype(int), 0, len(elrho)-1)
lp = LD + elpr[jT, iD]                      # log10 p = log10 rho + log10(p/rho)
# the EOS table is NaN where no gas solution exists (below ~71 K); those nodes
# must not be substituted, or the (T,p) lookup below indexes on garbage
lp_ok = np.isfinite(lp)
lp = np.where(lp_ok, lp, cP[0])

# the ck mean at those (log T, log p), clamped to its grid; mark where it was clamped
x = np.clip((LT - cT[0])/(cT[1]-cT[0]), 0, len(cT)-1.001)
y = np.clip((lp - cP[0])/(cP[1]-cP[0]), 0, len(cP)-1.001)
i0 = x.astype(int)
j0 = y.astype(int)
fx, fy = x - i0, y - j0
ck_here = ((1-fx)*(1-fy)*cK[i0, j0] + fx*(1-fy)*cK[i0+1, j0]
           + (1-fx)*fy*cK[i0, j0+1] + fx*fy*cK[i0+1, j0+1])
logR = LD - 3.0*LT + 18.0            # the coordinate the stellar sources are tabulated in
inrange = ((LT >= cT[0]) & (LT <= cT[-1]) & (lp >= cP[0]) & (lp <= cP[-1])
           & (logR >= -8.0) & (logR <= 1.0) & lp_ok & np.isfinite(ck_here))

# weight: 0 = stellar, 1 = ck
w = np.zeros_like(sK)
lo, hi, hi2 = np.log10(a.tlo), np.log10(a.thi), np.log10(a.thi2)
w = np.where(LT < lo, np.clip((LT - (lo - 0.3))/0.3, 0, 1), 1.0)
w = np.where(LT > hi, np.clip((hi2 - LT)/(hi2 - hi), 0, 1), w)
w = np.where(inrange, w, 0.0)

out = (1.0 - w)*sK + w*np.where(inrange, ck_here, sK)
use = inrange & (w > 0.5)
d = np.abs(sK - ck_here)[use]
print('where the ck mean is used, |log kappa_ck - log kappa_stellar| median %.3f, '
      '90%% %.3f, max %.3f over %d nodes' % (np.median(d), np.percentile(d, 90),
                                             d.max(), d.size))
for tl, th in [(3.3, 3.5), (3.5, 3.7), (3.7, 3.8)]:
    m = use & (LT >= tl) & (LT < th)
    if m.sum():
        print('   %5.0f-%5.0f K: median %+0.3f dex over %d nodes'
              % (10**tl, 10**th, np.median((ck_here - sK)[m]), m.sum()))
print('ck fraction: %.1f%% of nodes take it fully, %.1f%% partially'
      % (100*(w > 0.999).mean(), 100*((w > 0) & (w < 0.999)).mean()))
lt0, dlt, ld0, dld = grid
with open(a.out, 'w') as f:
    f.write('# AthenaK stellar Rosseland opacity table, log10 kappa_R [cm^2/g] on '
            '(log10 T[K], log10 rho[g/cm3])\n')
    f.write('# %s spliced with the CORRELATED-K derived mean (%s) over T = %g..%g K\n'
            '# and only where BOTH are valid (logR in [-8,1]); blended out to %g K\n'
            '# above. Grains keep the stellar table at the cold end.\n'
            % (a.stellar.split('/')[-1], a.ck.split('/')[-1], a.tlo, a.thi, a.thi2))
    f.write('# grid: nT nD lTmin dlT lDmin dlD\n')
    f.write('# %d %d %g %g %g %g\n' % (len(sT), len(sD), lt0, dlt, ld0, dld))
    f.write('# then nT*nD rows, T slowest: log10 kappa_R\n')
    for i in range(len(sT)):
        for j in range(len(sD)):
            f.write('%.5f\n' % out[i, j])
print('wrote', a.out, out.shape)
