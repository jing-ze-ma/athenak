"""ledger.py: close the energy/mass ledger of budget_rot300/ON (ebudget.txt, problem/budget_dt)
and compare with the hst booking (dhj.user.hst, dhj.hydro.hst).  usage: python3 ledger.py [ARM]"""
import sys
import numpy as np
D = sys.argv[1] if len(sys.argv) > 1 else 'ON'
PR = 1.101535e5
NQ = 20
N = ['FXE', 'FXM', 'HDVE', 'HDVM', 'SRCE', 'SRCM', 'FLRE', 'FLRM', 'QCK', 'FIR', 'QSW',
     'TOPST', 'TOPIO', 'FLRS', 'CKfloor', 'CKdemax', 'CKother', 'CKlincap', 'CKlin',
     'CNT', 'E', 'M']
blocks, cur, tt = [], None, []
for ln in open(D + '/ebudget.txt'):
    if ln.startswith('#'):
        if cur:
            blocks.append(np.array(cur))
        cur = []
        tt.append(float(ln.split()[2]))
    else:
        cur.append([float(v) for v in ln.split()[1:]])
blocks.append(np.array(cur))
A = np.array(blocks)          # (ndump, n=nx1+1, NQ+2)
t = np.array(tt)
a, b = 0, len(t) - 1
T = t[b] - t[a]
d = (A[b] - A[a]).T           # (NQ+2, n)
n1 = d.shape[1] - 1
pb = np.load('shellp.npy')
q = {k: d[i] for i, k in enumerate(N)}
dE = q['E'][:n1]
dM = q['M'][:n1]
fe = q['FXE']
fm = q['FXM']
pred = fe[:n1] - fe[1:] - q['HDVE'][:n1] + q['SRCE'][:n1] + q['FLRE'][:n1] + q['QCK'][:n1]
predm = fm[:n1] - fm[1:] - q['HDVM'][:n1] + q['SRCM'][:n1] + q['FLRM'][:n1]
res, resm = dE - pred, dM - predm
ckx = q['QCK'][:n1] - (q['FIR'][:n1] - q['FIR'][1:] + q['QSW'][:n1])
print('# %s: %d dumps, t %.6e..%.6e = %.3f rot' % (D, len(t), t[a], t[b], T/PR))
print('# GLOBAL ENERGY [erg/s] (time means over the window)')
g = dict(dE=dE.sum(), Fbot=fe[0], Ftop=-fe[n1], HDV=-q['HDVE'][:n1].sum(),
         SRC=q['SRCE'][:n1].sum(), FLR=q['FLRE'][:n1].sum(), FLRsplit=q['FLRS'][:n1].sum(),
         QCK=q['QCK'][:n1].sum(), RES=res.sum(), CKX=ckx.sum(),
         ck_fluxbooked=q['FIR'][0] - q['FIR'][n1] + q['QSW'][:n1].sum(),
         CKfloor=q['CKfloor'][:n1].sum(), CKdemax=q['CKdemax'][:n1].sum(),
         CKother=q['CKother'][:n1].sum(), CKlincap=q['CKlincap'][:n1].sum(),
         CKlin=q['CKlin'][:n1].sum(), FIRtop=q['FIR'][n1], FIRbot=q['FIR'][0], QSW=q['QSW'][:n1].sum())
for k, v in g.items():
    print('  %-14s %+.4e' % (k, v/T))
print('  counts (cells x dumps? cumulative) floor/demax/other/lincap/lin/full/lin:', q['CNT'][:7])
ts = q['TOPST']
io = q['TOPIO']
print('  top face: applied FXE %+.4e  stage1 %+.4e  stage2 %+.4e  (erg/s, outward +)'
      % (fe[n1]/T, ts[0]/T, ts[2]/T))
print('  top face: applied FXM %+.4e  stage1 %+.4e  stage2 %+.4e  (g/s)'
      % (fm[n1]/T, ts[1]/T, ts[3]/T))
print('  top inflow faces E %+.4e M %+.4e | outflow faces E %+.4e M %+.4e'
      % (io[0]/T, io[1]/T, io[2]/T, io[3]/T))
print('# GLOBAL MASS [g/s]: dM %+.4e Fbot %+.4e Ftop %+.4e HDV %+.4e SRC %+.4e FLR %+.4e RES %+.4e'
      % (dM.sum()/T, fm[0]/T, -fm[n1]/T, -q['HDVM'][:n1].sum()/T, q['SRCM'][:n1].sum()/T,
         q['FLRM'][:n1].sum()/T, resm.sum()/T))
# hst booking over the same window
u = np.loadtxt(D + '/dhj.user.hst')
h = np.loadtxt(D + '/dhj.hydro.hst')
lab = [s.split('=')[1] for s in open(D + '/dhj.user.hst').readlines()[1].split('[')[1:]]
lab = [s.split()[0] for s in lab]
s = (u[:, 0] > t[a] + 1) & (u[:, 0] <= t[b] + 1)
um = {k: u[s, i].mean() for i, k in enumerate(lab)}
s10 = s & (np.abs(((u[:, 0] - 3.304605e7)/(PR/10)) - np.round((u[:, 0] - 3.304605e7)/(PR/10))) < 1e-3)
print('# HST means (%d rows; %d at the production 10/rot cadence):' % (s.sum(), s10.sum()))
for k in ['Lir_top', 'Lsw_abs', 'Lrad_bot', 'Etot_top', 'Etot_bot', 'Mdot_top', 'Efloor',
          'Mfloor', 'Efloor_rt']:
    print('  %-10s %+.4e  (10/rot subset %+.4e)' % (k, um[k], u[s10, lab.index(k)].mean()))
fl_top = um['Etot_top'] - um['Lir_top'] + um['Lsw_abs']
print('  hst fluid top (Etot_top - Lir_top + Lsw_abs) %+.4e' % fl_top)
predh = um['Etot_bot'] - um['Etot_top'] + um['Efloor']
sh = (h[:, 0] >= t[a] - 1) & (h[:, 0] <= t[b] + 1)
c = np.polyfit(h[sh, 0], h[sh, 6], 1)
print('  hst-booked dE/dt %+.4e ; hydro.hst tot-E slope %+.4e ; endpoints %+.4e ; ledger dE %+.4e'
      % (predh, c[0], (h[sh, 6][-1] - h[sh, 6][0])/(h[sh, 0][-1] - h[sh, 0][0]), dE.sum()/T))
print('  gap (measured - hst-booked) %+.4e ; mass gap dM - (-Mdot_top + Mfloor) %+.4e'
      % (dE.sum()/T - predh, dM.sum()/T - (-um['Mdot_top'] + um['Mfloor'])))
# per-shell table for the top 20 shells and the group sums
print('# PER SHELL (erg/s): i p[bar] dE Fin(i)-Fout(i+1) -HDV SRC FLR QCK CKX RES | dM RESM')
for i in list(range(n1)):
    print('%3d %.2e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e | %+.3e %+.3e' % (
        i, pb[i], dE[i]/T, (fe[i] - fe[i+1])/T, -q['HDVE'][i]/T, q['SRCE'][i]/T,
        q['FLRE'][i]/T, q['QCK'][i]/T, ckx[i]/T, res[i]/T, dM[i]/T, resm[i]/T))
print('# GROUPS (erg/s): p range  dE  -HDV  SRC  FLR  QCK  CKX  RES | FXE at group top face')
edges = [1e9, 100, 10, 1, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0]
for hi, lo in zip(edges[:-1], edges[1:]):
    m = (pb < hi) & (pb >= lo)
    if not m.any():
        continue
    it = np.where(m)[0].max() + 1
    print('  %7.0e-%7.0e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e | %+.3e' % (
        lo, hi, dE[m].sum()/T, -q['HDVE'][:n1][m].sum()/T, q['SRCE'][:n1][m].sum()/T,
        q['FLRE'][:n1][m].sum()/T, q['QCK'][:n1][m].sum()/T, ckx[m].sum()/T, res[m].sum()/T,
        fe[it]/T))
# per-interval global gap, to see if steady
print('# per-dump-interval: t_rot  dE/dt  hst-style(FXE_bot-FXE_top+ckbooked+FLR)  QCK-ckbooked')
for k in range(1, len(t)):
    dd = (A[k] - A[k-1]).T
    dt = t[k] - t[k-1]
    ckb = dd[9][0] - dd[9][n1] + dd[10][:n1].sum()
    print('  %.3f %+.3e %+.3e %+.3e' % ((t[k] - 3.304605e7)/PR, dd[20][:n1].sum()/dt,
          (dd[0][0] - dd[0][n1] + ckb + dd[6][:n1].sum())/dt, (dd[8][:n1].sum() - ckb)/dt))
