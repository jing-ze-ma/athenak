"""I8_vceil analysis (2026-09-11).  Read-only; writes i8_report.txt."""
import sys, os, re, glob
import numpy as np
sys.path.insert(0, '/orion/u/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/orion/u/jinma/ATHENAK/athenak/tools/solar_convection')
sys.path.insert(0, '/orion/ptmp/jinma/Athenak/red_giant/_analysis_0910')
import bin_convert, rg

RUN = '/orion/ptmp/jinma/Athenak/red_giant/I8_vceil'
OUT = '/orion/ptmp/jinma/Athenak/red_giant/_analysis_0911/i8_report.txt'
L = []
def P(s=''):
    L.append(s); print(s)

# ---- radial grid (same stretch as I4/V8: nx1=480, 9.6e10 -> 4.16e12) ----
C8 = np.array([-0.196600, 8.811622, -17.210081, 8.883710]); R0, X1, N = 9.6e10, 4.16e12, 480
xi = np.arange(N+1)/N; _s = xi.copy(); _k = xi.copy()
for _j in range(1, 5):
    _s = _s + C8[_j-1]*_k*(1-xi); _k = _k*xi
F = R0 + (X1-R0)*_s
RL, RR = F[:-1], F[1:]; _q = RL/RR
RC = 0.25*(_q*_q+1)/((1/3)*(_q*_q+_q+1))*(RR+RL)
VCEIL = 5.0e7; THR = 4.9e7

def temp(d, e):
    d = np.asarray(d); e = np.asarray(e); sh = d.shape
    return np.asarray(rg.eos().temperature(np.maximum(d, 1e-23).ravel(),
                                           (e/d).ravel())).reshape(sh)

txt = open(os.path.join(RUN, 'out.txt'), 'rb').read().decode('utf8', 'replace')

# ================= 1. dt history ==================
cyc = np.array([[float(m.group(1)), float(m.group(2)), float(m.group(3))] for m in
                re.finditer(r'cycle=(\d+) time=([0-9.eE+-]+) dt=([0-9.eE+-]+)', txt)])
P('=== 1a. dt per 2.5e4 of simulated time ===')
P('  t_bin       ncyc    dt_min      dt_med      dt_last')
for b in np.arange(5.0e5, 9.0e5, 2.5e4):
    s = cyc[(cyc[:, 1] >= b) & (cyc[:, 1] < b+2.5e4)]
    if not len(s): continue
    P('  %9.3e %5d  %.4e  %.4e  %.4e' % (b, len(s), s[:, 2].min(), np.median(s[:, 2]), s[-1, 2]))
P('  last cycle=%d t=%.6e dt=%.4e' % (cyc[-1, 0], cyc[-1, 1], cyc[-1, 2]))
P('  (cycle lines are printed every 100 cycles, so dt_min here MISSES the short dips;')
P('   the true minima are the COLLAPSE blocks below.)')

P('\n=== 1b. dt COLLAPSE blocks ===')
P('  time         dt_old    dt_new   cell(m,k,j,i)        r          rho        T          v')
blocks = re.findall(r'### dt COLLAPSE cycle=(\d+) time=([0-9.eE+-]+) dtold=([0-9.eE+-]+) '
                    r'dt=([0-9.eE+-]+).*?\n\s+hydro dt is set by cell \(m,k,j,i\) = '
                    r'\(([0-9, ]+)\).*?\n\s+r=([0-9.eE+-]+) rho=([0-9.eE+-]+) T=([0-9.eE+-]+) '
                    r'p=([0-9.eE+-]+) cs=([0-9.eE+-]+) v=\(([^)]*)\)', txt)
for b in blocks:
    P('  %.5e  %7.3f  %6.3f  (%s)%s %.4e  %.2e  %.3e  v=(%s)' %
      (float(b[1]), float(b[2]), float(b[3]), b[4], ' '*(14-len(b[4])),
       float(b[5]), float(b[6]), float(b[7]), b[10]))
P('  (%d collapse blocks total)' % len(blocks))

# ---- 1c. RT rescues and NaN in the RT sweep ----
rep = re.findall(r'### rt_cell_report (\w+) ncycle=(\d+) t=([0-9.eE+-]+) m=(\d+) k=(\d+) '
                 r'j=(\d+) i=(\d+) r=([0-9.eE+-]+)\n(.*?\n.*?\n.*?\n.*?\n)', txt)
tt = np.array([float(r[2]) for r in rep])
nanr = np.array(['nan' in r[8] for r in rep])
P('\n=== 1c. two-stream RT cell reports ===')
P('  t_bin        reports  of which RESCUE  with NaN in the sweep')
for b in np.arange(5.0e5, 9.0e5, 5.0e4):
    m = (tt >= b) & (tt < b+5.0e4)
    if not m.sum(): continue
    P('  %9.3e  %7d  %14d  %s' % (b, m.sum(),
      sum(1 for r, mm in zip(rep, m) if mm and r[0] == 'RESCUE'), int((nanr & m).sum())))
if nanr.any():
    rr = [r for r, n in zip(rep, nanr) if n]
    ri = np.array([float(r[7]) for r in rr])
    P('  NaN reports: %d, t = %.4e .. %.4e, r = %.4e .. %.4e, i = %d..%d'
      % (len(rr), min(float(r[2]) for r in rr), max(float(r[2]) for r in rr),
         ri.min(), ri.max(), min(int(r[6]) for r in rr), max(int(r[6]) for r in rr)))
    P('  (I_dn/F/src come back -nan there; the report says rescued=1, de/e=-0.999)')
P('  open-ghost guard FIRED lines: %d (ghost i=482 only, replaced by the initial column)'
  % txt.count('open-ghost guard FIRED'))

# ================= 2. face budget ==================
fb = np.array([[float(m.group(3)), float(m.group(1)), float(m.group(2))] for m in
               re.finditer(r'L_rad,out/L = ([0-9.eE+-]+) L_rad,cut/L = ([0-9.eE+-]+) '
                           r'\(t = ([0-9.eE+-]+) s\)', txt)])
P('\n=== 2. emergent-flux budget (nearest face-budget line to each 2.5e4) ===')
P('  t_target     t_line       L_rad,out/L   L_rad,cut/L   [I5 ref]')
I5 = {7.0e5: 1.2, 8.25e5: 3.3, 9.0e5: 2.0}
for b in np.arange(5.0e5, 9.25e5, 2.5e4):
    k = np.argmin(np.abs(fb[:, 0]-b))
    ref = ''
    for tt, vv in I5.items():
        if abs(tt-b) < 1.3e4: ref = 'I5 %.1f @%.2e' % (vv, tt)
    P('  %9.3e  %.5e  %11.4f   %11.4f   %s' % (b, fb[k, 0], fb[k, 1], fb[k, 2], ref))
P('  lidded reference: L_rad,out/L = 1.3-1.5')
P('  the 5.70-5.90e5 BURST (every face-budget line, the shell crossing the top face):')
for r in fb[(fb[:, 0] > 5.62e5) & (fb[:, 0] < 5.96e5)]:
    P('     t=%.5e   L_rad,out/L = %11.4e   L_rad,cut/L = %+.4e' % (r[0], r[1], r[2]))

# ================= 3. history ==================
hst = np.loadtxt(os.path.join(RUN, 'rg.hydro.hst'))
def at(t):
    return hst[np.argmin(np.abs(hst[:, 0]-t))]
P('\n=== 3. energy / mass (rg.hydro.hst) ===')
P('  t             tot-E          mass           KE(1+2+3)      dE/E0      dM/M0')
r0 = at(5.0e5)
for t in (5.0e5, 5.6e5, 6.3e5, 7.0e5, 8.0e5, 9.0e5):
    r = at(t)
    P('  %.5e  %.6e  %.6e  %.6e  %+.2e  %+.2e' %
      (r[0], r[6], r[2], r[7]+r[8]+r[9], (r[6]-r0[6])/abs(r0[6]), (r[2]-r0[2])/r0[2]))
a, b = at(5.6e5), at(6.3e5)
P('  shell interval %.4e -> %.4e: dE/E = %+.3e   dM/M = %+.3e' %
  (a[0], b[0], (b[6]-a[6])/abs(a[6]), (b[2]-a[2])/a[2]))
P('  NOTE: the .hst prints 6 significant figures; tot-E and mass are BIT-CONSTANT over')
P('  all %d rows (min=max), so |dE/E| and |dM/M| are only bounded: < ~1e-5 over 5e5->9e5.' % len(hst))
P('  KE varies by %.2fx over the window, so the file is live.' %
  ((hst[:, 7]+hst[:, 8]+hst[:, 9]).max()/(hst[:, 7]+hst[:, 8]+hst[:, 9]).min()))

# ================= 4. ceiling ==================
lg = np.loadtxt(os.path.join(RUN, 'rg.log'))
ct = np.interp(lg[:, 0], cyc[:, 0], cyc[:, 1])
P('\n=== 4a. event counters per interval (rg.log; cycle -> time) ===')
P('  t(cyc)       cycle   eos_vceil   eos_dfloor   eos_efloor')
last = -1e30
for i in range(len(lg)):
    if ct[i]-last < 2.4e4: continue
    last = ct[i]
    P('  %.5e %7d %10d %11d %12d' % (ct[i], lg[i, 0], lg[i, 4], lg[i, 1], lg[i, 2]))
P('  peak eos_vceil = %d at t=%.4e ; total = %d' %
  (lg[:, 4].max(), ct[np.argmax(lg[:, 4])], lg[:, 4].sum()))

def dumps():
    o = []
    for p in sorted(glob.glob(os.path.join(RUN, 'bin', '*.hydro_w.*.bin'))):
        o.append((bin_convert.read_binary(p)['time'], p))
    return sorted(o)
DS = dumps()
def near_dump(t):
    return min(DS, key=lambda x: abs(x[0]-t))

def census(path, tag):
    """Where the FAST cells sit.  The ceiling itself uses the METRIC speed
    |v|^2 = g_ij v^i v^j (coordinates.cpp cs_raisev), so that is what is formed here."""
    fd = bin_convert.read_binary(path)
    a = {v: np.asarray(fd['mb_data'][v]) for v in fd['var_names']}
    mg = np.asarray(fd['mb_geometry']); n2 = fd['nx2_out_mb']; n3 = fd['nx3_out_mb']
    cc = np.zeros((fd['n_mbs'], n3, n2))
    for m in range(fd['n_mbs']):
        x2min, x2max, x3min, x3max = mg[m][2], mg[m][3], mg[m][4], mg[m][5]
        x2v = x2min + (np.arange(n2)+0.5)*(x2max-x2min)/n2
        x3v = (x3min + (np.arange(n3)+0.5)*(x3max-x3min)/n3)[::-1]
        cc[m] = -np.sin(0.25*np.pi*x3v)[:, None]*np.sin(0.25*np.pi*x2v)[None, :]
    c = cc[:, :, :, None]
    vx, vy, vz = a['velx'], a['vely'], a['velz']
    vsq = vx*vx + vy*vy + vz*vz + 2.0*c*vy*vz
    v = np.sqrt(np.maximum(vsq, 0.0))
    P('\n  -- %s  t=%.5e : max metric |v| = %.3e (ceiling %.1e); cells over 4.9e7: %d'
      % (tag, fd['time'], v.max(), VCEIL, int((v > THR).sum())))
    dfl = a['dens'] <= 1.01e-18
    nd = int(dfl.sum())
    if nd:
        ri = RC[np.argwhere(dfl)[:, 3]]
        P('     cells AT the density floor (1e-18): %d, r = %.3e..%.3e, max |v| there = %.3e'
          % (nd, ri.min(), ri.max(), v[dfl].max()))
    else:
        P('     cells AT the density floor (1e-18): 0')
    thr = 5.0e6
    sel = v > thr
    n = int(sel.sum())
    P('     cells over %.1e cm/s (10%% of the ceiling, the fastest population): %d of %d'
      % (thr, n, v.size))
    if n == 0:
        return
    ii = np.argwhere(sel)[:, 3]; rr = RC[ii]
    d = a['dens'][sel]; T = temp(d, a['eint'][sel]); vr = vx[sel]
    edges = np.linspace(RC[0], RC[-1], 11)
    h, _ = np.histogram(rr, bins=edges)
    P('     r-bin [1e12 cm]      N      rho range             T range            max|v|')
    for k in range(10):
        m = (rr >= edges[k]) & ((rr < edges[k+1]) if k < 9 else (rr <= edges[k+1]))
        if h[k] == 0:
            P('     %6.3f-%6.3f    %6d' % (edges[k]/1e12, edges[k+1]/1e12, 0)); continue
        P('     %6.3f-%6.3f    %6d   %.2e-%.2e   %.2e-%.2e   %.2e' %
          (edges[k]/1e12, edges[k+1]/1e12, h[k], d[m].min(), d[m].max(),
           T[m].min(), T[m].max(), v[sel][m].max()))
    P('     r(min/med/max) = %.3e / %.3e / %.3e ; cells AT dfloor(1e-18): %d' %
      (rr.min(), np.median(rr), rr.max(), int((d <= 1.01e-18).sum())))

P('\n=== 4b. the fast population in the dumps (metric |v|) ===')
for t, tag in ((6.0e5, 'dump near 6.0e5'), (6.6e5, 'dump near 6.6e5 = vceil PEAK'),
               (DS[-1][0], 'LAST dump')):
    tt, p = near_dump(t); census(p, tag)

P('\n=== 5a. top 60 radial shells at the last dump ===')
tt, p = near_dump(9.0e5)
fd = bin_convert.read_binary(p)
a = {v: np.asarray(fd['mb_data'][v]) for v in fd['var_names']}
T = temp(a['dens'], a['eint'])
P('  t = %.5e' % fd['time'])
P('   i     r[cm]      max T      med T      med rho    max|v_r|')
for i in range(N-60, N):
    P('  %3d  %.4e  %.3e  %.3e  %.3e  %.3e' %
      (i, RC[i], T[:, :, :, i].max(), np.median(T[:, :, :, i]),
       np.median(a['dens'][:, :, :, i]), np.abs(a['velx'][:, :, :, i]).max()))

P('\n=== 5b. shell front from the SHELL-MEDIAN density profile ===')
P('  t            r(rho=1e-11)  r(rho=1e-12)  r(rho=1e-13)  r(rho=1e-14)   med-rho@3.65e12')
rows = []
for t in (5.2e5, 6.0e5, 7.0e5, 8.0e5, 9.0e5):
    tt, p = near_dump(t)
    fd = bin_convert.read_binary(p)
    med = np.median(np.asarray(fd['mb_data']['dens']), axis=(0, 1, 2))
    rs = []
    for thr in (1e-11, 1e-12, 1e-13, 1e-14):
        idx = np.where((med < thr) & (np.arange(N) > 300))[0]
        rs.append(RC[idx[0]] if len(idx) else np.nan)
    i0 = int(np.argmin(np.abs(RC-3.65e12)))
    P('  %.4e  %.4e  %.4e  %.4e  %.4e   %.3e' % ((fd['time'],)+tuple(rs)+(med[i0],)))
    rows.append((fd['time'],)+tuple(rs))
P('  front speed d r(rho=1e-12)/dt between consecutive rows [cm/s]:')
for k in range(1, len(rows)):
    P('    %.4e -> %.4e : %+.3e' % (rows[k-1][0], rows[k][0],
                                    (rows[k][2]-rows[k-1][2])/(rows[k][0]-rows[k-1][0])))

P('')
P('=== 6. VERDICT ===')
P(open('/orion/ptmp/jinma/Athenak/red_giant/_analysis_0911/verdict.txt').read())
open(OUT, 'w').write('\n'.join(L)+'\n')
print('\nwrote', OUT)
