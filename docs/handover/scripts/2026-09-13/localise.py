#!/usr/bin/env python3
"""Part A: localise the sp MHD-module excess (mhd_b0 vs hyd_ctl) in the dumps."""
import sys, os
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc

ROOT = '/viper/u2/jinma/ATHENAK/bench/sp_excess'
OUT = os.path.join(ROOT, 'analysis')
NX1, NX2, NX3 = 128, 64, 128
R0, R1 = 9.44e9, 2.0556e10
CPOLY = [-0.068392, -2.191487, 2.464818, -1.366698]
ATHETA = 3.0
GRAV, AP = 942.0, 9.44e9        # g(r) = grav*(ap/r)^2 (point mass)
GM = GRAV*AP*AP

def stretch_r(x):
    xi = (x-R0)/(R1-R0); u = xi.copy(); xik = xi.copy()
    for k in range(1, 5):
        u += CPOLY[k-1]*xik*(1.0-xi); xik = xik*xi
    return R0 + (R1-R0)*u

def stretch_th(t):
    xi = t/np.pi
    return np.pi/2.0*(1.0+np.sinh(ATHETA*(2.0*xi-1.0))/np.sinh(ATHETA))

# --- global grid (exactly as Coordinates builds it) ---
x1f = stretch_r(np.linspace(R0, R1, NX1+1))
x2f = stretch_th(np.linspace(0.0, np.pi, NX2+1))
x3f = np.linspace(0.0, 2.0*np.pi, NX3+1)
x1v = 0.5*(x1f[:-1]+x1f[1:]); x2v = 0.5*(x2f[:-1]+x2f[1:]); x3v = 0.5*(x3f[:-1]+x3f[1:])
# volume element dV = (r_r^3-r_l^3)/3 * (cos th_l - cos th_r) * dphi
dV1 = (x1f[1:]**3 - x1f[:-1]**3)/3.0
dV2 = np.cos(x2f[:-1]) - np.cos(x2f[1:])
dV3 = np.diff(x3f)
dV = dV3[:, None, None]*dV2[None, :, None]*dV1[None, None, :]   # (k,j,i)

def load(arm, var, snap):
    f = bc.read_binary(f'{ROOT}/{arm}/out/bin/dhj.{var}.{snap:05d}.bin')
    n1, n2, n3 = f['nx1_mb'], f['nx2_mb'], f['nx3_mb']
    out = {}
    for v in ['dens', 'velx', 'vely', 'velz', 'eint']:
        g = np.zeros((NX3, NX2, NX1))
        d = np.asarray(f['mb_data'][v])
        for m, (l1, l2, l3, lev) in enumerate(f['mb_logical']):
            g[l3*n3:(l3+1)*n3, l2*n2:(l2+1)*n2, l1*n1:(l1+1)*n1] = d[m]
        out[v] = g
    out['time'] = f['time']
    return out

def derived(d):
    r = {}
    r['mass'] = d['dens']*dV
    r['eint'] = d['eint']*dV
    v2h = d['vely']**2 + d['velz']**2
    r['ke_h'] = 0.5*d['dens']*v2h*dV
    r['ke_r'] = 0.5*d['dens']*d['velx']**2*dV
    r['pe'] = -GM/x1v[None, None, :]*d['dens']*dV
    r['etot'] = r['eint'] + r['ke_h'] + r['ke_r'] + r['pe']
    return r

fh = open(os.path.join(OUT, 'localise.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); fh.write(s+'\n')

cosphi = np.cos(x3v)
dayk = cosphi > 0.0

P('='*86)
P('PART A -- localising the sp MHD-module excess:  mhd_b0  MINUS  hyd_ctl')
P('grid 128(r) x 64(theta) x 128(phi); r stretched poly, theta sinh a=3, phi uniform')
P('substellar at phi=0 (mu = sin(theta) cos(phi)), dayside = cos(phi)>0')
P('quantities: eint, KE_h = 0.5 rho (v_th^2+v_ph^2), KE_r, PE = -GM rho/r, etot = sum')
P('='*86)

snaps = [(1, '0.25 rot'), (2, '0.50 rot')]
store = {}
for snap, lbl in snaps:
    b = derived(load('mhd_b0', 'mhd_w_bcc', snap))
    h = derived(load('hyd_ctl', 'hydro_w', snap))
    store[snap] = (b, h)
    P('')
    P('#'*86)
    P(f'### {lbl}')
    P('#'*86)
    P('')
    P('--- GLOBAL (volume integrals) ---')
    P(f'{"qty":8s} {"mhd_b0":>14s} {"hyd_ctl":>14s} {"diff":>14s} {"ratio":>8s}')
    for q in ['mass', 'eint', 'ke_r', 'ke_h', 'pe', 'etot']:
        B, H = b[q].sum(), h[q].sum()
        P(f'{q:8s} {B:14.6e} {H:14.6e} {B-H:14.4e} {B/H if H!=0 else 0:8.3f}')

    # ---- (1) vs radius ----
    P('')
    P('--- (1) SHELL SUMS vs RADIUS INDEX i (diff = mhd_b0 - hyd_ctl) ---')
    de = (b['etot']-h['etot']).sum(axis=(0, 1))
    dk = (b['ke_h']-h['ke_h']).sum(axis=(0, 1))
    dm = (b['mass']-h['mass']).sum(axis=(0, 1))
    di = (b['eint']-h['eint']).sum(axis=(0, 1))
    dp = (b['pe']-h['pe']).sum(axis=(0, 1))
    kb = b['ke_h'].sum(axis=(0, 1)); kh = h['ke_h'].sum(axis=(0, 1))
    bands = [('deep   i 0-57', 0, 58), ('photo  i 58-98', 58, 99),
             ('top    i 99-127', 99, 128)]
    P(f'{"band":16s} {"r range [cm]":>24s} {"d_etot":>12s} {"d_eint":>12s} '
      f'{"d_PE":>12s} {"d_KEh":>12s} {"KEh ratio":>10s} {"d_mass":>12s}')
    for nm, a, z in bands:
        P(f'{nm:16s} {x1v[a]:11.4e}-{x1v[z-1]:11.4e} {de[a:z].sum():12.4e} '
          f'{di[a:z].sum():12.4e} {dp[a:z].sum():12.4e} {dk[a:z].sum():12.4e} '
          f'{kb[a:z].sum()/max(kh[a:z].sum(),1e-99):10.3f} {dm[a:z].sum():12.4e}')
    P('')
    P('  per-8-cell radial bins:')
    P(f'  {"i":>8s} {"r[cm]":>11s} {"d_etot":>12s} {"d_eint":>12s} {"d_KEh":>12s} '
      f'{"KEh ratio":>10s} {"d_mass":>12s}')
    for a in range(0, NX1, 8):
        z = a+8
        P(f'  {str(a)+"-"+str(z-1):>8s} {x1v[a]:11.4e} {de[a:z].sum():12.4e} '
          f'{di[a:z].sum():12.4e} {dk[a:z].sum():12.4e} '
          f'{kb[a:z].sum()/max(kh[a:z].sum(),1e-99):10.3f} {dm[a:z].sum():12.4e}')
    iw_e = int(np.argmax(np.abs(de))); iw_k = int(np.argmax(np.abs(dk)))
    iw_m = int(np.argmax(np.abs(dm)))
    P(f'  worst radius: etot i={iw_e} (r={x1v[iw_e]:.4e}, {de[iw_e]:.3e}); '
      f'KEh i={iw_k} ({dk[iw_k]:.3e}); mass i={iw_m} ({dm[iw_m]:.3e})')

    # ---- (2) vs latitude within band ----
    P('')
    P('--- (2) vs LATITUDE (theta row j), and day/night ---')
    for nm, a, z in bands:
        P(f'  [{nm}]')
        dej = (b['etot']-h['etot'])[:, :, a:z].sum(axis=(0, 2))
        dkj = (b['ke_h']-h['ke_h'])[:, :, a:z].sum(axis=(0, 2))
        dmj = (b['mass']-h['mass'])[:, :, a:z].sum(axis=(0, 2))
        kbj = b['ke_h'][:, :, a:z].sum(axis=(0, 2))
        khj = h['ke_h'][:, :, a:z].sum(axis=(0, 2))
        grp = [('j 0-1  (N pole)', 0, 2), ('j 2-3', 2, 4), ('j 4-7', 4, 8),
               ('j 8-15', 8, 16), ('j 16-31 (N mid/eq)', 16, 32),
               ('j 32-47 (S mid/eq)', 32, 48), ('j 48-55', 48, 56),
               ('j 56-59', 56, 60), ('j 60-61', 60, 62),
               ('j 62-63 (S pole)', 62, 64)]
        P(f'    {"rows":20s} {"theta[deg]":>12s} {"d_etot":>12s} {"d_KEh":>12s} '
          f'{"KEh ratio":>10s} {"d_mass":>12s}')
        for gn, ja, jz in grp:
            P(f'    {gn:20s} {np.degrees(x2v[ja]):5.2f}-{np.degrees(x2v[jz-1]):5.2f} '
              f'{dej[ja:jz].sum():12.4e} {dkj[ja:jz].sum():12.4e} '
              f'{kbj[ja:jz].sum()/max(khj[ja:jz].sum(),1e-99):10.3f} '
              f'{dmj[ja:jz].sum():12.4e}')
        dd = (b['etot']-h['etot'])[:, :, a:z]
        dk2 = (b['ke_h']-h['ke_h'])[:, :, a:z]
        P(f'    day  (cosphi>0): d_etot {dd[dayk].sum():12.4e}  '
          f'd_KEh {dk2[dayk].sum():12.4e}')
        P(f'    night(cosphi<0): d_etot {dd[~dayk].sum():12.4e}  '
          f'd_KEh {dk2[~dayk].sum():12.4e}')

    # ---- (3) polar-row anomaly at the worst radius ----
    P('')
    P(f'--- (3) POLAR-ROW / SEAM anomaly at worst-KEh radius i={iw_k} '
      f'(r={x1v[iw_k]:.4e}) ---')
    bd = load('mhd_b0', 'mhd_w_bcc', snap); hd = load('hyd_ctl', 'hydro_w', snap)
    def rowstat(d, i):
        vphi = np.abs(d['velz'][:, :, i]).mean(axis=0)
        vth = np.abs(d['vely'][:, :, i]).mean(axis=0)
        T = (d['eint'][:, :, i]/d['dens'][:, :, i]).mean(axis=0)   # ~ propto T
        rho = d['dens'][:, :, i].mean(axis=0)
        return vphi, vth, T, rho
    for i in sorted(set([iw_k, iw_e, iw_m])):
        P(f'  radius i={i} r={x1v[i]:.4e}')
        vb = rowstat(bd, i); vh = rowstat(hd, i)
        P(f'    {"j":>4s} {"theta":>7s} {"|vphi|_b":>11s} {"|vphi|_h":>11s} {"rat":>6s}'
          f' {"|vth|_b":>11s} {"|vth|_h":>11s} {"rat":>6s} {"e/rho_b":>11s}'
          f' {"e/rho_h":>11s} {"rat":>6s} {"rho ratio":>9s}')
        for j in list(range(0, 5))+list(range(28, 36))+list(range(59, 64)):
            P(f'    {j:4d} {np.degrees(x2v[j]):7.2f} {vb[0][j]:11.4e} {vh[0][j]:11.4e}'
              f' {vb[0][j]/max(vh[0][j],1e-99):6.2f} {vb[1][j]:11.4e} {vh[1][j]:11.4e}'
              f' {vb[1][j]/max(vh[1][j],1e-99):6.2f} {vb[2][j]:11.4e} {vb[2][j]*0+vh[2][j]:11.4e}'
              f' {vb[2][j]/max(vh[2][j],1e-99):6.2f} {vb[3][j]/max(vh[3][j],1e-99):9.4f}')
        # pole-adjacent vs next-2 rows
        for pn, p2, n2 in [('N', [0, 1], [2, 3]), ('S', [62, 63], [60, 61])]:
            rv = vb[0][p2].mean()/max(vh[0][p2].mean(), 1e-99)
            rn = vb[0][n2].mean()/max(vh[0][n2].mean(), 1e-99)
            tv = vb[2][p2].mean()/max(vh[2][p2].mean(), 1e-99)
            tn = vb[2][n2].mean()/max(vh[2][n2].mean(), 1e-99)
            P(f'    {pn} pole: |vphi| ratio(b/h) pole-2rows {rv:.3f} vs next-2 {rn:.3f}'
              f' ; e/rho ratio {tv:.4f} vs {tn:.4f}')
        # phi-seam check: is the diff concentrated at any k?
        dk3 = np.abs(bd['velz'][:, :, i]-hd['velz'][:, :, i]).mean(axis=1)
        km = int(np.argmax(dk3))
        P(f'    phi structure of |dvphi| (mean over j): max at k={km} '
          f'(phi={np.degrees(x3v[km]):.1f} deg) = {dk3[km]:.4e}, '
          f'median = {np.median(dk3):.4e}, ratio {dk3[km]/max(np.median(dk3),1e-99):.2f}')
        P(f'    k=0 {dk3[0]:.4e}  k=1 {dk3[1]:.4e}  k=127 {dk3[127]:.4e}  '
          f'k=126 {dk3[126]:.4e}  (a phi seam would spike at k=0/127)')

    # ---- (4) where does the mass go ----
    P('')
    P('--- (4) MASS GAIN: where ---')
    dmf = b['mass']-h['mass']
    P(f'  total d_mass {dmf.sum():.4e}')
    P(f'  inner 4 radial rows i=0-3   : {dmf[:, :, 0:4].sum():.4e}')
    P(f'  i=4-57                      : {dmf[:, :, 4:58].sum():.4e}')
    P(f'  photosphere i=58-98         : {dmf[:, :, 58:99].sum():.4e}')
    P(f'  top i=99-127                : {dmf[:, :, 99:128].sum():.4e}')
    P(f'  outer 4 rows i=124-127      : {dmf[:, :, 124:128].sum():.4e}')
    P(f'  polar rows j=0,1,62,63 (all i): {dmf[:, [0,1,62,63], :].sum():.4e}')
    P(f'  j=2..61                       : {dmf[:, 2:62, :].sum():.4e}')
    dmj2 = dmf.sum(axis=(0, 2))
    P('  d_mass per theta row (worst 8):')
    for j in np.argsort(-np.abs(dmj2))[:8]:
        P(f'    j={j:3d} theta={np.degrees(x2v[j]):6.2f}  {dmj2[j]:12.4e}')

P('')
P('='*86)
fh.close()

# ---- maps ----
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b, h = store[2]
    bd = load('mhd_b0', 'mhd_w_bcc', 2); hd = load('hyd_ctl', 'hydro_w', 2)
    de = (b['etot']-h['etot'])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    m0 = de.sum(axis=0).T          # (i,j)
    v = np.nanpercentile(np.abs(m0), 99)
    p = ax[0].pcolormesh(np.degrees(x2v), np.arange(NX1), m0, cmap='RdBu_r',
                         vmin=-v, vmax=v)
    ax[0].set_xlabel('theta [deg]'); ax[0].set_ylabel('radial index i')
    ax[0].set_title('d(etot) mhd_b0-hyd_ctl, phi-summed, 0.5 rot')
    fig.colorbar(p, ax=ax[0])
    m1 = (b['ke_h']-h['ke_h']).sum(axis=0).T
    v = np.nanpercentile(np.abs(m1), 99)
    p = ax[1].pcolormesh(np.degrees(x2v), np.arange(NX1), m1, cmap='RdBu_r',
                         vmin=-v, vmax=v)
    ax[1].set_xlabel('theta [deg]'); ax[1].set_ylabel('radial index i')
    ax[1].set_title('d(KE_h), phi-summed, 0.5 rot')
    fig.colorbar(p, ax=ax[1])
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'localise_rtheta.png'), dpi=110)

    iw = 110
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.0))
    for n, (dd, t) in enumerate([(bd['velz'][:, :, iw].T, 'mhd_b0 v_phi'),
                                 (hd['velz'][:, :, iw].T, 'hyd_ctl v_phi'),
                                 ((bd['velz']-hd['velz'])[:, :, iw].T, 'diff v_phi')]):
        v = np.nanpercentile(np.abs(dd), 99)
        p = ax[n].pcolormesh(np.degrees(x3v), np.degrees(x2v), dd, cmap='RdBu_r',
                             vmin=-v, vmax=v)
        ax[n].set_title(f'{t}  i={iw}'); ax[n].set_xlabel('phi [deg]')
        ax[n].set_ylabel('theta [deg]'); ax[n].invert_yaxis()
        fig.colorbar(p, ax=ax[n])
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'localise_vphi_map.png'), dpi=110)
    print('maps written')
except Exception as e:
    print('plot skipped:', e)
