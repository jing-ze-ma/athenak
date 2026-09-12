"""RG_fofc_long deathmap: geography of the HOT cells and of the FOFC flags in the
ghost-inclusive hydro_w / hydro_fofc dumps of the deterministically reproduced dt
collapse (slurm 11649850).  Analysis only: writes only into this directory.

Grid: cubed sphere, 6 panels x 32x32 columns, nx1 = 480 radial (stretched),
24 MeshBlocks of 480x16x16, nghost = 3.  With ghost_zones=true a block array is
(k,j,i) = (22,22,486); the ACTIVE range is 3..18 in k,j and 3..482 in i, i.e. the
array index equals the CODE index that the dt-collapse message prints.
j <-> x2 = xi, k <-> x3 = eta, panel = gid//4.
"""
import sys, os, glob
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert

B    = os.environ.get('DEATHMAP_B', '/viper/u2/jinma/ATHENAK/bench/RG_fofc_long/deathmap/')
OUT  = os.environ.get('DEATHMAP_OUT', B + 'analysis/')
NP, NC, NR, NG, NB = 6, 32, 480, 3, 16
R0, R1 = 9.6e10, 4.16e12
CPOLY  = np.array([-0.196600, +8.811622, -17.210081, +8.883710])
GAM, MU = 5.0/3.0, 0.60
MH, KB  = 1.6726219e-24, 1.380649e-16
TFAC = (MU*MH/KB)*(GAM-1.0)                 # T = TFAC * e/rho
DEATH = (15, 9, 18, 339)                    # gid,k,j,i  CODE indices from out.txt
RSTAR, RJOIN = 3.2e12, 3.652820e12

# ---------------- radial grid (same map as coordinates.cpp cs_coord1d_1) --------------
def stretch(r, c=CPOLY, r0=R0, r1=R1):
    xi = (np.asarray(r, float)-r0)/(r1-r0)
    u = xi.copy(); xik = xi.copy()
    for k in range(1, len(c)+1):
        u = u + c[k-1]*xik*(1.0-xi); xik = xik*xi
    return r0 + (r1-r0)*u
XE  = stretch(np.linspace(R0, R1, NR+1))
Q   = XE[:-1]/XE[1:]
RVa = 0.25*(Q*Q+1.0)/((1.0/3.0)*(Q*Q+Q+1.0))*(XE[1:]+XE[:-1])
DRa = XE[1:]-XE[:-1]
RV  = np.empty(NR+2*NG); RV[NG:NG+NR] = RVa
for g in range(NG):
    RV[NG-1-g]   = RVa[0]  - (g+1)*DRa[0]
    RV[NG+NR+g]  = RVa[-1] + (g+1)*DRa[-1]

def load(path):
    r = bin_convert.read_binary(path)
    return (r['time'], r['cycle'],
            {v: np.asarray(r['mb_data'][v]).astype(np.float64) for v in r['var_names']},
            np.array(r['mb_logical']), list(r['var_names']))

def geom(nk, nj, mbl):
    ng = (nk-NB)//2
    kk = np.arange(nk)[:, None]-ng
    jj = np.arange(nj)[None, :]-ng
    J = np.empty((24, nk, nj), int); K = np.empty((24, nk, nj), int)
    for m in range(24):
        J[m] = NB*mbl[m, 1]+jj        # xi  index (j axis)
        K[m] = NB*mbl[m, 2]+kk        # eta index (k axis)
    AG = (J < 0) | (J >= NC) | (K < 0) | (K >= NC)
    Jc, Kc = np.clip(J, 0, NC-1), np.clip(K, 0, NC-1)
    DE = np.minimum(np.minimum(Jc, NC-1-Jc), np.minimum(Kc, NC-1-Kc))
    DV = np.maximum(np.minimum(Jc, NC-1-Jc), np.minimum(Kc, NC-1-Kc))
    return J, K, AG, DE, DV

RBINS = [(0, 3.0e12), (3.0e12, 3.3e12), (3.3e12, 3.42e12), (3.42e12, 3.55e12),
         (3.55e12, RJOIN), (RJOIN, 3.9e12), (3.9e12, 1e13)]
RBNAME = ['<3.0e12', '3.0-3.3', '3.3-3.42', '3.42-3.55', '3.55-join', 'join-3.9', '>3.9e12']

def rbin(r):
    out = np.zeros(r.shape, int)
    for q, (a, b) in enumerate(RBINS):
        out[(r >= a) & (r < b)] = q
    return out

def SUB(tser):
    o = tser[::max(1, len(tser)//12)]
    if o[-1] is not tser[-1]: o = o + [tser[-1]]
    return o

def main():
    log = open(OUT+'deathmap.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); log.write(s+'\n'); log.flush()

    wf = sorted(glob.glob(B+'bin/rg.hydro_w.*.bin'))
    ff = sorted(glob.glob(B+'bin/rg.hydro_fofc.*.bin'))
    P('RG_fofc_long DEATHMAP  --  %d hydro_w dumps, %d hydro_fofc dumps' % (len(wf), len(ff)))
    P('T_proxy = (mu m_H/k)(gamma-1) e/rho,  mu=0.60, gamma=5/3  ->  T = %.4e * (e/rho)' % TFAC)
    P('(the tabulated EOS is not available here; T_proxy is the ionised-H ideal-gas value.')
    P(' In the cool envelope it OVERSTATES T by the ionisation/H2 heat capacity; at the')
    P(' runaway cell -- fully ionised, 1e8 K+ -- it is accurate to tens of percent.)')
    P('')

    t0, c0, d0, mbl, vn = load(wf[0])
    nk, nj, ni = d0['dens'].shape[1:]
    ngA, ngR = (nk-NB)//2, (ni-NR)//2
    J, K, AG, DE, DV = geom(nk, nj, mbl)
    RVi = RV if ngR == NG else RVa
    rgh = np.zeros(ni, bool)
    if ngR: rgh[:ngR] = True; rgh[ni-ngR:] = True
    P('hydro_w block shape (k,j,i) = (%d,%d,%d); angular ghosts %d, radial ghosts %d; vars %s'
      % (nk, nj, ni, ngA, ngR, vn))
    gd, kd, jd, idd = DEATH        # CODE indices (is=js=ks=3)
    # prefer the cell the run itself printed (out.txt of THIS reproduction)
    import re
    try:
        txt = open(B+'out.txt', errors='replace').read()
        mm = re.findall(r'hydro dt is set by cell \(m,k,j,i\) = \((\d+),(\d+),(\d+),(\d+)\)'
                        r' gid = (\d+)', txt)
        if mm:
            a = mm[0]
            kd, jd, idd, gd = int(a[1]), int(a[2]), int(a[3]), int(a[4])
            P('out.txt: %d dt-collapse cell reports.  ALL OF THEM, with geography:' % len(mm))
            P('  %4s %4s %3s %3s %4s %11s %4s %4s %8s' %
              ('n','gid','k','j','i','r[cm]','dE','dV','panel'))
            seen = []
            for a2 in mm:
                key = (int(a2[4]), int(a2[1]), int(a2[2]), int(a2[3]))
                if key not in [x[0] for x in seen]: seen.append([key, 0])
                for x in seen:
                    if x[0] == key: x[1] += 1
            for key, n in seen:
                g2, k2, j2, i2 = key
                ka, ja, ia = k2-(NG-ngA), j2-(NG-ngA), i2-(NG-ngR)
                P('  %4d %4d %3d %3d %4d %11.5e %4d %4d %8d' %
                  (n, g2, k2, j2, i2, RVi[ia], DE[g2,ka,ja], DV[g2,ka,ja], g2//4))
    except OSError:
        pass
    kd -= (NG-ngA); jd -= (NG-ngA); idd -= (NG-ngR)   # -> array indices of this dump
    P('DEATH CELL (m,k,j,i)=(3,%d,%d,%d) rank 1 -> gid %d = panel %d, block (lx2,lx3)=(%d,%d)'
      % (kd, jd, idd, gd, gd//4, mbl[gd,1], mbl[gd,2]))
    P('   global column (J=xi,K=eta) = (%d,%d);  dist to panel edge %d cells, to cube vertex'
      ' %d cells;  r = %.5e cm = %.4f R*' %
      (J[gd,kd,jd], K[gd,kd,jd], DE[gd,kd,jd], DV[gd,kd,jd], RVi[idd], RVi[idd]/RSTAR))
    P('')

    def b4(a): return np.broadcast_to(a[:, :, :, None], (24, nk, nj, ni))
    AG4, DE4, DV4 = b4(AG), b4(DE), b4(DV)
    RG4 = np.broadcast_to(rgh[None, None, None, :], (24, nk, nj, ni))
    ACT = ~(AG4 | RG4)
    R4  = np.broadcast_to(RVi[None, None, None, :], (24, nk, nj, ni))
    RB4 = rbin(R4)
    EB4 = np.minimum(DE4, 3)
    # reference population of ACTIVE cells per edge bin / radial bin
    refE = np.array([int(((EB4 == q) & ACT).sum()) for q in range(4)], float)
    refR = np.array([int(((RB4 == q) & ACT).sum()) for q in range(len(RBINS))], float)
    NACT = refE.sum()
    P('ACTIVE cells %d;  by distance from a panel edge: dE=0 %d (%.2f%%), 1 %d, 2 %d, 3+ %d'
      % (NACT, refE[0], 100*refE[0]/NACT, refE[1], refE[2], refE[3]))
    P('')

    # storage for the time series
    tser = []
    S = {}      # dump tag -> dict of scalars

    # ============ pass over the hydro_w dumps ============
    P('=== 1. HOTTEST ACTIVE CELL PER DUMP  (whole domain, and restricted to r > 3.3e12) ===')
    P('%5s %11s %7s | %10s %10s %4s %3s %3s %4s %10s %2s %2s | %10s %10s %2s %2s | %7s %7s'
      % ('dump','time','cycle','max e/rho','T[K]','gid','k','j','i','r','dE','dV',
         'T out[K]','r out','dE','dV','n>1e8','n>1e12'))
    rows1, rows1b, rows4 = [], [], []
    outer = ACT & (R4 > 3.3e12)
    for f in wf:
        t, c, d, _, _ = load(f)
        er = d['eint']/d['dens']; T = TFAC*er
        Ta = np.where(ACT, T, -np.inf); ix = np.unravel_index(Ta.argmax(), T.shape)
        To = np.where(outer, T, -np.inf); ox = np.unravel_index(To.argmax(), T.shape)
        n8  = int(((T > 1e8) & ACT).sum()); n12 = int(((T > 1e12) & ACT).sum())
        P('%5s %11.5e %7d | %10.3e %10.3e %4d %3d %3d %4d %10.4e %2d %2d | %10.3e %10.4e %2d %2d | %7d %7d'
          % (f[-9:-4], t, c, er[ix], T[ix], ix[0], ix[1], ix[2], ix[3], RVi[ix[3]],
             DE[ix[0],ix[1],ix[2]], DV[ix[0],ix[1],ix[2]],
             T[ox], RVi[ox[3]], DE[ox[0],ox[1],ox[2]], DV[ox[0],ox[1],ox[2]], n8, n12))
        # geography rows
        r1 = {}
        for thr, nm in ((1e8,'1e8'), (1e12,'1e12')):
            m = (T > thr) & ACT; n = int(m.sum())
            r1[nm] = (n, [int(((EB4 == q) & m).sum()) for q in range(4)],
                      [int(((DV4 <= 1) & m).sum()), int(((DV4 <= 2) & m).sum())],
                      [int(((RB4 == q) & m).sum()) for q in range(len(RBINS))],
                      (R4[m].min(), R4[m].max()) if n else (0, 0))
        rows1b.append((f[-9:-4], t, r1))
        # ghost census
        m8 = T > 1e8
        rows4.append((f[-9:-4], t, int((m8 & ACT).sum()), int((m8 & AG4 & ~RG4).sum()),
                      int((m8 & RG4 & ~AG4).sum()), int((m8 & AG4 & RG4).sum()),
                      float(np.where(ACT, T, -np.inf).max()),
                      float(np.where(~ACT, T, -np.inf).max())))
        # time series at the death cell + 3x3 angular neighbourhood at i = idd
        sl = (gd, slice(kd-1, kd+2), slice(jd-1, jd+2), idd)
        tser.append(dict(tag=f[-9:-4], t=t, c=c,
                         er=float(er[gd,kd,jd,idd]), rho=float(d['dens'][gd,kd,jd,idd]),
                         vr=float(d['velx'][gd,kd,jd,idd]),
                         vh=float(np.hypot(d['vely'][gd,kd,jd,idd], d['velz'][gd,kd,jd,idd])),
                         ei=float(d['eint'][gd,kd,jd,idd]),
                         er9=er[sl].copy(), rho9=d['dens'][sl].copy(),
                         vr9=d['velx'][sl].copy(),
                         vh9=np.hypot(d['vely'][sl], d['velz'][sl]),
                         # radial neighbourhood along the column
                         ercol=er[gd, kd, jd, idd-3:idd+4].copy(),
                         rhocol=d['dens'][gd, kd, jd, idd-3:idd+4].copy(),
                         vrcol=d['velx'][gd, kd, jd, idd-3:idd+4].copy()))
        del d, er, T
    P('')

    P('=== 1b. GEOGRAPHY of the hot ACTIVE cells: distance in cells from the nearest panel edge ===')
    P('    (a dE=0 cell is ON a panel seam; dV<=1 means within 1 cell of a CUBE VERTEX)')
    for nm, thr in (('1e8', 1e8), ('1e12', 1e12)):
        P('-- T_proxy > %s K --   expected shares if uniform: dE=0 %.3f  1 %.3f  2 %.3f  3+ %.3f'
          % (nm, *(refE/NACT)))
        P('%5s %11s %8s | %6s %6s %6s %6s | %6s %6s | %6s | %10s %10s' %
          ('dump','time','n','dE=0','dE=1','dE=2','dE>=3','dV<=1','dV<=2','enr0','r min','r max'))
        for tag, t, r1 in rows1b:
            n, eb, dv, rb, rr = r1[nm]
            if n == 0:
                P('%5s %11.5e %8d | %6s' % (tag, t, 0, '-')); continue
            enr0 = (eb[0]/n)/(refE[0]/NACT)
            P('%5s %11.5e %8d | %6d %6d %6d %6d | %6d %6d | %6.2f | %10.4e %10.4e'
              % (tag, t, n, eb[0], eb[1], eb[2], eb[3], dv[0], dv[1], enr0, rr[0], rr[1]))
        P('   radial histogram of the same cells (bins %s):' % ' | '.join(RBNAME))
        for tag, t, r1 in rows1b:
            n, eb, dv, rb, rr = r1[nm]
            if n: P('   %5s %11.5e  %s' % (tag, t, ' '.join('%7d' % x for x in rb)))
        P('')

    P('=== 4. ARE THE HOT CELLS IN THE GHOST LAYERS? (T_proxy > 1e8 K, all cells) ===')
    P('%5s %11s | %9s %9s %9s %6s | %11s %11s' %
      ('dump','time','active','ang ghost','rad ghost','both','max T act','max T ghost'))
    for row in rows4:
        P('%5s %11.5e | %9d %9d %9d %6d | %11.4e %11.4e' % row)
    P('')

    # ============ 3. time series at the death cell ============
    P('=== 3. TIME SERIES at the death cell gid %d (k,j,i)=(%d,%d,%d), r=%.4e ==='
      % (gd, kd, jd, idd, RVi[idd]))
    P('%5s %11s %7s | %11s %11s %11s | %11s %11s | %11s' %
      ('dump','time','cycle','rho','e/rho','T_proxy[K]','v_r','|v_h|','eint'))
    for s in tser:
        P('%5s %11.5e %7d | %11.4e %11.4e %11.4e | %11.3e %11.3e | %11.4e' %
          (s['tag'], s['t'], s['c'], s['rho'], s['er'], TFAC*s['er'], s['vr'], s['vh'], s['ei']))
    P('')
    P('e-folding rates over each dump interval (dt = t_n - t_{n-1}):')
    P('%5s %11s | %9s %9s %9s | %11s %11s %11s' %
      ('dump','time','dln e/dt','dln rho/dt','dln(e/rho)','tau_e[s]','tau_rho[s]','tau_er[s]'))
    for a, b in zip(tser[:-1], tser[1:]):
        dt = b['t']-a['t']
        def rate(x, y):
            return (np.log(max(y,1e-300))-np.log(max(x,1e-300)))/dt
        re_, rr_, rer = rate(a['ei'],b['ei']), rate(a['rho'],b['rho']), rate(a['er'],b['er'])
        P('%5s %11.5e | %9.2e %9.2e %9.2e | %11.4e %11.4e %11.4e' %
          (b['tag'], b['t'], re_, rr_, rer,
           1/re_ if re_ else np.inf, 1/rr_ if rr_ else np.inf, 1/rer if rer else np.inf))
    P('')
    P('3x3 ANGULAR neighbourhood (k-1..k+1 rows, j-1..j+1 cols; j=%d is the LAST ACTIVE'
      ' xi row = the panel seam, j=%d is the FIRST ANGULAR GHOST) at i=%d:' % (jd, jd+1, idd))
    for s in SUB(tser):
        P('  dump %s  t=%.5e' % (s['tag'], s['t']))
        for nm, a in (('T_proxy', TFAC*s['er9']), ('rho', s['rho9']), ('v_r', s['vr9'])):
            P('    %-8s %s' % (nm, ' | '.join(' '.join('%10.3e' % v for v in row) for row in a)))
    P('')
    P('RADIAL neighbourhood i=%d..%d of the death column:' % (idd-3, idd+3))
    for s in SUB(tser):
        P('  dump %s t=%.5e  T %s' % (s['tag'], s['t'],
          ' '.join('%9.3e' % v for v in TFAC*s['ercol'])))
        P('                          rho %s' % ' '.join('%9.3e' % v for v in s['rhocol']))
        P('                          v_r %s' % ' '.join('%9.3e' % v for v in s['vrcol']))
    P('')
    np.save(OUT+'tseries.npy', np.array(tser, dtype=object), allow_pickle=True)

    # ============ 2. FOFC map ============
    if ff:
        t, c, d, mbl2, vn2 = load(ff[0])
        fv = vn2[0]
        P('=== 2. FOFC MAP  (variable %r; per-cell flag count since the previous dump) ===' % fv)
        nk2, nj2, ni2 = d[fv].shape[1:]
        P('hydro_fofc block shape (k,j,i) = (%d,%d,%d)' % (nk2, nj2, ni2))
        P('%5s %11s %7s | %9s %9s %9s | %6s %6s %6s %6s | %6s %6s | %10s %4s %3s %3s %4s %10s' %
          ('dump','time','cycle','n flagged','sum flags','max flags','dE=0','dE=1','dE=2','dE>=3',
           'dV<=1','dV<=2','enr(dE=0)','gid','k','j','i','r'))
        frows = []
        for f in ff:
            t, c, d, _, _ = load(f)
            a = d[fv]
            m = (a > 0) & ACT
            n = int(m.sum())
            ix = np.unravel_index(np.where(ACT, a, -1).argmax(), a.shape)
            if n == 0:
                P('%5s %11.5e %7d | %9d %9s' % (f[-9:-4], t, c, 0, '-')); frows.append((f[-9:-4],t,None)); continue
            eb = [int(((EB4 == q) & m).sum()) for q in range(4)]
            enr = (eb[0]/n)/(refE[0]/NACT)
            P('%5s %11.5e %7d | %9d %9.0f %9.0f | %6d %6d %6d %6d | %6d %6d | %10.2f %4d %3d %3d %4d %10.4e'
              % (f[-9:-4], t, c, n, float(a[ACT].sum()), float(a[ACT].max()),
                 eb[0], eb[1], eb[2], eb[3],
                 int(((DV4 <= 1) & m).sum()), int(((DV4 <= 2) & m).sum()), enr,
                 ix[0], ix[1], ix[2], ix[3], RVi[ix[3]]))
            frows.append((f[-9:-4], t, dict(mask_rb=[int(((RB4 == q) & m).sum()) for q in range(len(RBINS))])))
            del d, a, m
        P('   radial histogram of the flagged cells (bins %s):' % ' | '.join(RBNAME))
        for tag, t, r in frows:
            if r: P('   %5s %11.5e  %s' % (tag, t, ' '.join('%7d' % x for x in r['mask_rb'])))
        P('')

        # coincidence hot <-> flagged, matched by dump tag
        P('=== 2b. DO THE FLAGGED CELLS COINCIDE WITH THE HOT CELLS? ===')
        P('%5s %11s | %8s %8s %8s | %8s %8s | %10s' %
          ('dump','time','n hot','n flag','n both','P(flag|hot)','P(hot|flag)','flag@death'))
        wtimes = {round(ss['t'], 3): fw for ss, fw in zip(tser, wf)}
        for f in ff:
            tag = f[-9:-4]
            t, c, d, _, _ = load(f); a = d[fv]
            key = min(wtimes, key=lambda x: abs(x-t))
            if abs(key-t) > 1.0:
                del d, a; continue
            _, _, dw, _, _ = load(wtimes[key])
            T = TFAC*dw['eint']/dw['dens']
            hot = (T > 1e8) & ACT; flg = (a > 0) & ACT
            nb = int((hot & flg).sum()); nh = int(hot.sum()); nf = int(flg.sum())
            P('%5s %11.5e | %8d %8d %8d | %8.4f %8.4f | %10.0f' %
              (tag, t, nh, nf, nb, nb/nh if nh else 0.0, nb/nf if nf else 0.0,
               float(a[gd, kd, jd, idd])))
            del d, dw, a, T, hot, flg
    log.close()

if __name__ == '__main__':
    main()
