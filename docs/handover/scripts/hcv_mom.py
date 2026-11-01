"""hclean_1031 HCV: deep horizontal 1-2-1 filter of the momentum (rho v) on a remap.dat
(RADREMPE) of an identity C32 remap.  usage: hcv_mom.py IN.dat OUT.dat X1F_DUMP [CFULL CZERO]

- per radial cell i: Cartesian momentum per cell P = rho v dV filtered with the global
  shell smoother M_c (hodge_filter.Shell: tensor 1-2-1 across block and panel seams,
  doubly stochastic -> sum P conserved exactly), weight w(i) = 1 (i <= CFULL), cos^2 taper
  to 0 at CZERO (same taper as the B filter).
- then per shell one solid-rotation increment dOm about z plus one uniform velocity U,
  dP_c = rho_c dV_c (dOm z x r_c + U), solved so that sum P and L_z are both restored
  exactly (curved-sphere averaging changes L_z at O((dx/R)^2) otherwise).
- rho untouched; e_new = e + KE_old - KE_new per cell (KE with the code's cubed-sphere
  metric, kinetic() of dhj_remap), so E_tot = e + KE + rho Phi is conserved cell by cell.
"""
import os
import sys
import struct
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import hodge_filter as hf  # noqa: E402
import dhj_remap as dr  # noqa: E402
import bin_convert as bc  # noqa: E402

fin, fout, fdump = sys.argv[1:4]
cfull, czero = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (12, 20)
b = open(fin, 'rb').read()
assert b[:8] == b'RADREMPE'
nmb, n1, n2, n3 = struct.unpack('<4i', b[8:24])
rl = 20 + 5*n1*n2*n3*8
lls, recs = [], []
for m in range(nmb):
    o = 24 + m*rl
    lls.append(np.frombuffer(b[o:o + 20], '<i4').copy())
    recs.append(np.frombuffer(b[o + 20:o + rl], '<f8').reshape(5, n3, n2, n1).copy())
nb = 1 + max(ll[1] for ll in lls)
N = n2*nb
# panel arrays (6, N, N, n1)
Q = np.zeros((5, 6, N, N, n1))
for ll, rc in zip(lls, recs):
    _, lx2, lx3, _, pn = ll
    Q[:, pn, lx3*n3:(lx3 + 1)*n3, lx2*n2:(lx2 + 1)*n2] = rc
rho, u1, u2, u3, ei = Q
e = np.asarray(bc.read_binary(fdump)['mb_x1f'])[0].astype(np.float64)
assert len(e) == n1 + 1
ro = dr.centroids(e)
om, _, _ = dr.cs_shell_geom(N)
dV = om[None, :, :, None]*((e[1:]**3 - e[:-1]**3)/3.0)[None, None, None, :]
c = dr.cs_cosc(N)[None, :, :, None]
# frames
cq = 0.25*np.pi*(-1.0 + 2.0*(np.arange(N) + 0.5)/N)
XI, ET = np.meshgrid(cq, cq)
hh = 1e-7
X = np.stack([hf.pmap(p, XI, ET) for p in range(6)])
ex = np.stack([hf.pmap(p, XI + hh, ET) - hf.pmap(p, XI - hh, ET) for p in range(6)])
et = np.stack([hf.pmap(p, XI, ET + hh) - hf.pmap(p, XI, ET - hh) for p in range(6)])
ex /= np.linalg.norm(ex, axis=-1)[..., None]
et /= np.linalg.norm(et, axis=-1)[..., None]
cchk = np.abs((ex*et).sum(-1) - c[..., 0]).max()
print('frame check |e_xi.e_eta - cos_cell| max %.2e' % cchk)
det = 1.0 - c*c
w2 = (u2 - c*u3)/det           # contravariant
w3 = (u3 - c*u2)/det
V = (u1[..., None]*X[:, :, :, None, :] + w2[..., None]*ex[:, :, :, None, :]
     + w3[..., None]*et[:, :, :, None, :])          # (6, N, N, n1, 3) Cartesian v
KE0 = dr.kinetic(rho, rho*u1, rho*u2, rho*u3, c)
print('KE check |0.5 rho |v|^2 - kinetic()| / max %.2e'
      % (np.abs(0.5*rho*(V**2).sum(-1) - KE0).max()/KE0.max()))
sh = hf.Shell(N)
wc = hf.taper(n1, cfull, czero)
act = np.where(wc > 0)[0]
P = (rho*dV)[..., None]*V                            # (6, N, N, n1, 3)
Pc = P.reshape(sh.nc, n1, 3)
mdv = (rho*dV).reshape(sh.nc, n1)
Rc = X.reshape(sh.nc, 3)
zr = np.stack([-Rc[:, 1], Rc[:, 0], 0*Rc[:, 0]], -1)     # z x rhat
Pn = Pc.copy()
rows = []
for i in act:
    q = Pc[:, i]
    qn = q + wc[i]*(sh.Mc @ q - q)
    r = ro[i]
    lz = lambda pp: (r*(Rc[:, 0]*pp[:, 1] - Rc[:, 1]*pp[:, 0])).sum()  # noqa: E731
    dlz, dP = lz(q) - lz(qn), q.sum(0) - qn.sum(0)
    # unknowns dOm, Ux, Uy, Uz: dP_c = m_c (dOm r zr_c + U)
    bz = (mdv[:, i][:, None]*r*zr)                       # dP per unit dOm
    A = np.zeros((4, 4)); rhs = np.zeros(4)
    A[0, 0] = lz(bz); rhs[0] = dlz
    for k in range(3):
        e_k = np.zeros(3); e_k[k] = 1
        bu = mdv[:, i][:, None]*e_k[None, :]
        A[0, 1 + k] = lz(bu)
        A[1 + k, 0] = bz[:, k].sum()
        A[1 + k, 1 + k] = mdv[:, i].sum()
        rhs[1 + k] = dP[k]
    sol = np.linalg.solve(A, rhs)
    qn = qn + sol[0]*bz + mdv[:, i][:, None]*sol[1:][None, :]
    Pn[:, i] = qn
    rows.append((i, wc[i], (lz(qn) - lz(q))/max(abs(lz(q)), 1e-300),
                 np.abs(qn.sum(0) - q.sum(0)).max()/np.abs(q).sum(0).max(),
                 sol[0], np.abs(sol[1:]).max()))
Vn = (Pn/mdv[..., None]).reshape(V.shape)
u1n = (Vn*X[:, :, :, None, :]).sum(-1)
u2n = (Vn*ex[:, :, :, None, :]).sum(-1)
u3n = (Vn*et[:, :, :, None, :]).sum(-1)
keep = wc == 0      # cells outside the filter: bitwise copies (no frame round trip)
u1n[..., keep], u2n[..., keep], u3n[..., keep] = u1[..., keep], u2[..., keep], u3[..., keep]
Vn[..., keep, :] = V[..., keep, :]
KE1 = dr.kinetic(rho, rho*u1n, rho*u2n, rho*u3n, c)
KE1[..., keep] = KE0[..., keep]
ein = ei + (KE0 - KE1)
print('eint min %.4e -> %.4e; min (e_new / e_old) %.6f' % (ei.min(), ein.min(),
                                                        (ein/ei).min()))
print('\n i  w     dLz/Lz     |dP|/sum|P|   dOm[1/s]   |U|max[cm/s]')
for r_ in rows:
    print('%2d %.3f %+.2e %.2e %+.3e %.3e' % r_)
# totals
Lz = lambda VV: ((rho*dV)*ro[None, None, None, :]*(X[..., 0][..., None]*VV[..., 1]  # noqa
                                                    - X[..., 1][..., None]*VV[..., 0])).sum()
Pt = lambda VV: ((rho*dV)[..., None]*VV).sum(axis=(0, 1, 2, 3))  # noqa: E731
print('\ntotal L_z %.12e -> %.12e (rel %+.2e)' % (Lz(V), Lz(Vn), Lz(Vn)/Lz(V) - 1))
print('total P', Pt(V), '->', Pt(Vn), ' |dP|/sum|P| %.2e'
      % (np.abs(Pt(Vn) - Pt(V)).max()/np.abs(P).sum()))
E0 = ((ei + KE0)*dV).sum(); E1 = ((ein + KE1)*dV).sum()
print('sum (e + KE) dV %.12e -> %.12e (rel %+.2e); cellwise max |d(e+KE)|/(e+KE) %.2e'
      % (E0, E1, E1/E0 - 1, (np.abs((ein + KE1) - (ei + KE0))/(ei + KE0)).max()))
# KE removed and v noise per group
lat = np.degrees(np.arcsin(X[..., 2]))
eq = np.abs(lat) < 30
print('\ngroup  KE eq old [erg]  removed  | KE all old  removed')
for lo, hi in ((0, 1), (2, 10), (11, 20), (21, 35), (36, n1 - 1)):
    sl = slice(lo, hi + 1)
    k0, k1 = (KE0*dV)[..., sl], (KE1*dV)[..., sl]
    print('%3d-%-3d %.4e %+.4f | %.4e %+.4f' % (lo, hi, k0[eq].sum(), k1[eq].sum()/k0[eq].sum()
                                             - 1, k0.sum(), k1.sum()/k0.sum() - 1))
# write
Qn = np.stack([rho, u1n, u2n, u3n, ein])
with open(fout, 'wb') as fo:
    fo.write(b[:24])
    for ll, rc in zip(lls, recs):
        _, lx2, lx3, _, pn = ll
        fo.write(ll.astype('<i4').tobytes())
        fo.write(np.ascontiguousarray(Qn[:, pn, lx3*n3:(lx3 + 1)*n3,
                                         lx2*n2:(lx2 + 1)*n2], '<f8').tobytes())
# untouched-above check
d = max(np.abs(Qn[k][..., czero:] - Q[k][..., czero:]).max() for k in range(5))
print('max |change| in cells >= %d: %.2e' % (czero, d))
print('wrote', fout)
