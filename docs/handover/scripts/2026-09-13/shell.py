"""RG_fofc: shell-launch / vertex-chimney analysis.  Read-only; writes analysis/shell.txt."""
import sys, os, glob
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert

RUN = '/viper/u2/jinma/ATHENAK/bench/RG_fofc'
OUT = os.path.join(RUN, 'analysis', 'shell.txt')
L = []
def P(s=''):
    L.append(s); print(s)

# radial grid (same stretch as the input)
C8 = np.array([-0.196600, 8.811622, -17.210081, 8.883710]); R0, X1, N = 9.6e10, 4.16e12, 480
xi = np.arange(N+1)/N; _s = xi.copy(); _k = xi.copy()
for _j in range(1, 5):
    _s = _s + C8[_j-1]*_k*(1-xi); _k = _k*xi
F = R0 + (X1-R0)*_s
RL, RR = F[:-1], F[1:]
_q = RL/RR
RC = 0.25*(_q*_q+1)/((1/3)*(_q*_q+_q+1))*(RR+RL)
DVR = (RR**3 - RL**3)/3.0          # r^2 dr integrated exactly

JOIN = 3.6528e12
IJ = int(np.argmin(np.abs(RC-JOIN)))

# ---- solid angle per angular cell (equiangular gnomonic, x in [-1,1] -> a = pi/4 x) ----
def domega(x2l, x2r, x3l, x3r, n=8):
    """exact-ish dOmega = int dX dY /(1+X^2+Y^2)^{3/2}, X=tan(pi/4 x2)."""
    g, w = np.polynomial.legendre.leggauss(n)
    a2 = 0.5*(x2r-x2l)*g + 0.5*(x2l+x2r); w2 = 0.5*(x2r-x2l)*w
    a3 = 0.5*(x3r-x3l)*g + 0.5*(x3l+x3r); w3 = 0.5*(x3l+x3r)*0 + 0.5*(x3r-x3l)*w
    X = np.tan(0.25*np.pi*a2); Y = np.tan(0.25*np.pi*a3)
    JX = 0.25*np.pi/np.cos(0.25*np.pi*a2)**2
    JY = 0.25*np.pi/np.cos(0.25*np.pi*a3)**2
    f = (1.0 + X[None, :]**2 + Y[:, None]**2)**-1.5
    return float((w3*JY)[:, None].T @ np.zeros(0)) if False else float(
        np.sum((w3*JY)[:, None]*(w2*JX)[None, :]*f))

def dumps():
    o = []
    for p in sorted(glob.glob(os.path.join(RUN, 'bin', '*.hydro_w.*.bin'))):
        o.append((bin_convert.read_binary(p)['time'], p))
    return sorted(o)
DS = dumps()
def near(t):
    return min(DS, key=lambda x: abs(x[0]-t))

fd0 = bin_convert.read_binary(DS[0][1])
NMB = fd0['n_mbs']; N2 = int(fd0['nx2_out_mb']); N3 = int(fd0['nx3_out_mb'])
mg = np.asarray(fd0['mb_geometry']); lg = np.asarray(fd0['mb_logical'])

# per-cell solid angle, and the (k,j) index of the cube-vertex column of each block.
# array axes are (m, k, j, i); the k (x3) axis is REVERSED w.r.t. mb_geometry
# (see red-giant-vertex-chimney: m=0 vertex cell is (k,j)=(nx3-1,0) for the (0,0) block).
DOM = np.zeros((NMB, N3, N2))
for m in range(NMB):
    x2l, x2r, x3l, x3r = mg[m][2], mg[m][3], mg[m][4], mg[m][5]
    e2 = np.linspace(x2l, x2r, N2+1)
    e3 = np.linspace(x3l, x3r, N3+1)[::-1]     # reversed k
    for k in range(N3):
        for j in range(N2):
            DOM[m, k, j] = domega(e2[j], e2[j+1], min(e3[k], e3[k+1]), max(e3[k], e3[k+1]))
P('# solid-angle check: sum dOmega = %.8f  (4 pi = %.8f)' % (DOM.sum(), 4*np.pi))

VERT = []          # (m,k,j) of the 24 cube-vertex columns
for m in range(NMB):
    l2, l3 = lg[m][1], lg[m][2]
    j = 0 if l2 == 0 else N2-1
    k = N3-1 if l3 == 0 else 0
    VERT.append((m, k, j))
vmask = np.zeros((NMB, N3, N2), bool)
for m, k, j in VERT:
    vmask[m, k, j] = True
P('# %d cube-vertex columns identified (one per MeshBlock, 6 panels x 4 corners)' % vmask.sum())
P('# join r = %.5e  -> i = %d (RC[i] = %.5e)' % (JOIN, IJ, RC[IJ]))
P('')

VOL = DOM[:, :, :, None]*DVR[None, None, None, :]

# ================== A. mass above the join vs t ==================
P('=== A. mass above the join r=3.6528e12 (i>=%d) vs t ===' % IJ)
P('  t            M(>join)[g]   M(>3.8e12)    M(>4.0e12)    M_tot[g]     f(>join)')
I38 = int(np.argmin(np.abs(RC-3.8e12))); I40 = int(np.argmin(np.abs(RC-4.0e12)))
hist = []
for t, p in DS:
    fd = bin_convert.read_binary(p)
    d = np.asarray(fd['mb_data']['dens'])
    dm = d*VOL
    mt = dm.sum(); mj = dm[:, :, :, IJ:].sum()
    hist.append((t, mj, dm[:, :, :, I38:].sum(), dm[:, :, :, I40:].sum(), mt))
for r in hist:
    P('  %.4e  %.5e  %.5e  %.5e  %.5e  %.4e' % (r[0], r[1], r[2], r[3], r[4], r[1]/r[4]))
P('  (M_tot from the dump should equal the .hst mass 5.04768e34 to the accuracy of the')
P('   solid-angle quadrature; the RATIO column and the time trend are what matter.)')
P('')

# ================== B. mean v_r profile near the top ==================
P('=== B. shell-mean and shell-median v_r near the top ===')
for t in (5e5, 6e5, 7e5, 8e5, 9e5):
    tt, p = near(t)
    fd = bin_convert.read_binary(p)
    d = np.asarray(fd['mb_data']['dens']); vx = np.asarray(fd['mb_data']['velx'])
    w = DOM[:, :, :, None]
    P('  -- t = %.5e' % tt)
    P('     i     r[cm]      <v_r>_area   <v_r>_mass   med v_r     max v_r    min v_r    med rho')
    for i in list(range(300, 480, 20))+[470, 475, 479]:
        va = (vx[:, :, :, i]*DOM).sum()/DOM.sum()
        vm = (vx[:, :, :, i]*DOM*d[:, :, :, i]).sum()/((DOM*d[:, :, :, i]).sum())
        P('     %3d  %.4e  %+.3e   %+.3e   %+.3e  %+.3e  %+.3e  %.3e' %
          (i, RC[i], va, vm, np.median(vx[:, :, :, i]), vx[:, :, :, i].max(),
           vx[:, :, :, i].min(), np.median(d[:, :, :, i])))
P('')

# ================== C. the 24 cube-vertex columns vs the bulk ==================
P('=== C. cube-vertex columns vs bulk: R_rho = median(rho_vertex)/median(rho_bulk) ===')
P('  (chimney radii per red-giant-vertex-chimney: i ~ 295-320)')
P('  t           R_rho@295  R_rho@300  R_rho@310  R_rho@320  R_rho@396  minR(295-320)  vr_vtx@310  vr_bulk@310')
for t, p in DS:
    fd = bin_convert.read_binary(p)
    d = np.asarray(fd['mb_data']['dens']); vx = np.asarray(fd['mb_data']['velx'])
    rr = []
    for i in range(290, 330):
        rv = np.median(d[:, :, :, i][vmask]); rb = np.median(d[:, :, :, i][~vmask])
        rr.append(rv/rb)
    rr = np.array(rr)
    g = lambda i: rr[i-290]
    P('  %.4e  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f@i=%d  %+.3e  %+.3e' %
      (t, g(295), g(300), g(310), g(320),
       np.median(d[:, :, :, 396][vmask])/np.median(d[:, :, :, 396][~vmask]),
       rr[5:31].min(), 295+int(np.argmin(rr[5:31])),
       np.median(vx[:, :, :, 310][vmask]), np.median(vx[:, :, :, 310][~vmask])))
P('')
P('  -- last dump: the 24 vertex columns individually at i=310 --')
fd = bin_convert.read_binary(DS[-1][1])
d = np.asarray(fd['mb_data']['dens']); vx = np.asarray(fd['mb_data']['velx'])
rb = np.median(d[:, :, :, 310][~vmask])
P('     m   k   j   rho          rho/rho_bulk   v_r')
for m, k, j in VERT:
    P('     %2d %3d %3d  %.4e   %8.4f     %+.3e' % (m, k, j, d[m, k, j, 310],
                                                    d[m, k, j, 310]/rb, vx[m, k, j, 310]))
P('')

# ================== D. history ==================
P('=== D. rg.hydro.hst: KE, total E, mass ===')
hst = np.loadtxt(os.path.join(RUN, 'rg.hydro.hst'))
ke = hst[:, 7]+hst[:, 8]+hst[:, 9]
P('  t             tot-E          mass           KE_tot         KE1          KE2          KE3')
for r, k in zip(hst, ke):
    if abs(r[0] % 5e4) < 1.1e4 or r[0] == hst[-1, 0]:
        P('  %.5e  %.6e  %.6e  %.4e  %.4e  %.4e  %.4e'
          % (r[0], r[6], r[2], k, r[7], r[8], r[9]))
P('  tot-E  min/max = %.6e / %.6e  (spread %.2e)' % (hst[:, 6].min(), hst[:, 6].max(),
                                                     (hst[:, 6].max()-hst[:, 6].min())/hst[:, 6].mean()))
P('  mass   min/max = %.6e / %.6e  (spread %.2e)' % (hst[:, 2].min(), hst[:, 2].max(),
                                                     (hst[:, 2].max()-hst[:, 2].min())/hst[:, 2].mean()))
P('  KE_tot min/max = %.4e / %.4e ; last 1e5: %.4e -> %.4e' %
  (ke.min(), ke.max(), ke[hst[:, 0] >= 8.0e5][0], ke[-1]))
open(OUT, 'w').write('\n'.join(L)+'\n')
print('\nwrote', OUT)
