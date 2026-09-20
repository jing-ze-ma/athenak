"""tests_r24: energy budget of tests_r22/lr1 in 0.25-turnover windows. Read-only."""
import re, struct, sys
import numpy as np

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
A = B + '/tests_r22/lr1'
TURN = 4705.0
RS = 2.3717e11
OMEGA = 0.17678
LW = 4.0776e37
LSTAR = LW/OMEGA
TAUD = 500.0
TDAMP_END = 18820.0

# ---------- hst ----------
h = np.loadtxt(A + '/he4.hydro.hst')
th, dth, mass, mom1 = h[:, 0], h[:, 1], h[:, 2], h[:, 3]
etot = h[:, 6]
ke = h[:, 7] + h[:, 8] + h[:, 9]

# ---------- log face budget ----------
pat = re.compile(r'face budget: gain/L\s+inner=(\S+) outer=(\S+) \| erg in=(\S+) out=(\S+)'
                 r' \| g in=(\S+) out=(\S+)(?: \| L_rad,out/L = (\S+) L_rad,cut/L = (\S+))?'
                 r' \(t = (\S+) s\)')
rows = []
for line in open(B + '/tests_r22/lr1.log', errors='ignore'):
    m = pat.search(line)
    if m:
        g = m.groups()
        rows.append([float(g[8]), float(g[2]), float(g[3]), float(g[4]), float(g[5]),
                     float(g[6]) if g[6] else np.nan, float(g[7]) if g[7] else np.nan])
F = np.array(rows)
tf, Ein, Eout, Min, Mout, Lro, Lrc = F.T
print("log: %d budget lines, t %.1f -> %.1f s (%.3f -> %.3f turn)"
      % (len(tf), tf[0], tf[-1], tf[0]/TURN, tf[-1]/TURN))

# ---------- rt_profile ----------
def rd(fn):
    out = []
    f = open(fn, 'rb')
    while True:
        hh = f.read(16)
        if len(hh) < 16:
            break
        t, n1, nv = struct.unpack('<dii', hh)
        b = f.read(8*(n1 + nv*n1))
        if len(b) < 8*(n1 + nv*n1):
            break
        a = np.frombuffer(b, '<f8')
        out.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1).copy()))
    return out

P = rd(A + '/rt_profile.bin')
r = P[0][1]
rf = np.empty(len(r)+1)
rf[1:-1] = 0.5*(r[1:] + r[:-1])
rf[0] = r[0] - (rf[1]-r[0])
rf[-1] = r[-1] + (r[-1]-rf[-2])
vol = OMEGA*4*np.pi/3*(rf[1:]**3 - rf[:-1]**3)
tp = np.array([x[0] for x in P])
print("rt_profile: %d records, t %.1f -> %.1f (%.3f -> %.3f turn), n1=%d nv=%d"
      % (len(P), tp[0], tp[-1], tp[0]/TURN, tp[-1]/TURN, len(r), P[0][2].shape[0]))

# damping sink power at each profile time: sum_shell m_shell <v1>^2 / tau
# (the operator does m1 -> m1 - (1-exp(-dt/tau)) rho <v1>, and subtracts the KE change
#  from IEN; to leading order the power removed is sum rho V <v1>^2/tau ... but the
#  cell KE change is 0.5*((m1-g*rho<v>)^2 - m1^2)/rho, averaged over the shell this is
#  -(g)<v1>*<rho v1>/... -> evaluate exactly with the shell mean only, per-cell
#  fluctuations cancel at first order.)
Pdamp = np.array([ (P[i][2][0]*vol*P[i][2][1]**2).sum()/TAUD for i in range(len(P)) ])
Pdamp[tp >= TDAMP_END] = 0.0

def interp(t, x, y):
    return np.interp(t, x, y)

# ---------- windows ----------
print()
print("all powers in units of L_w = %.4e erg/s   (L_star = %.4e)" % (LW, LSTAR))
print()
hdr = ("%7s %7s | %9s %9s | %8s %8s %8s | %8s %8s | %8s | %8s | %9s"
       % ("t0", "t1", "E0", "E1", "dE/dt", "L_in", "-L_out", "hyd_in", "hyd_out",
          "damp", "resid", "dM/dt[g/s]"))
print(hdr)
print("turn    turn    | erg       erg       | L_w      L_w      L_w      | L_w"
      "      L_w      | L_w      | L_w      |")
tend = min(th[-1], tf[-1])
edges = np.arange(2.0, tend/TURN + 1e-9, 0.25)
if edges[-1] < tend/TURN - 0.02:
    edges = np.append(edges, tend/TURN)
res = []
for a, b in zip(edges[:-1], edges[1:]):
    t0, t1 = a*TURN, b*TURN
    dt = t1 - t0
    E0 = interp(t0, th, etot); E1 = interp(t1, th, etot)
    dEdt = (E1-E0)/dt
    M0 = interp(t0, th, mass); M1 = interp(t1, th, mass)
    dMdt = (M1-M0)/dt
    lout = np.nanmean(Lro[(tf >= t0) & (tf < t1)])*LSTAR
    hin = (interp(t1, tf, Ein) - interp(t0, tf, Ein))/dt
    hout = (interp(t1, tf, Eout) - interp(t0, tf, Eout))/dt
    sel = (tp >= t0) & (tp < t1)
    pd = Pdamp[sel].mean() if sel.any() else 0.0
    resid = dEdt - (LW - lout + hin + hout - pd)
    res.append((a, b, E0, E1, dEdt, LW, -lout, hin, hout, -pd, resid, dMdt))
    print("%7.2f %7.2f | %9.3e %9.3e | %8.3f %8.3f %8.3f | %8.3f %8.3f | %8.3f |"
          " %8.3f | %9.2e"
          % (a, b, E0, E1, dEdt/LW, 1.0, -lout/LW, hin/LW, hout/LW, -pd/LW,
             resid/LW, dMdt))
np.save('/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r24/budget.npy', np.array(res))
