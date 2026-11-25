"""p23.py STAR a-b: v_r and comoving T_eff at tau_R = 2/3 of each column, page definition (pages_c1_1124/src/reduce.py: tau at cell centres from the top,
linear in tau between the two cells around 2/3, quantity linear between cell centres), no smoothing. Columns without tau >= 2/3: NaN. read only."""
import sys, os, types, re, time
import numpy as np
sys.modules.setdefault('h5py', types.ModuleType('h5py')); sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d'); import bin_convert as bc
G = '/resnick/groups/carnegie_poc/jingze/'; T2 = G + 'fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
R = dict(bsg=('prod_c1_1124/bsg/', 'bsg3d', 0.62, G + 'fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt'), m35mms=('prod_c1_1124/m35mms/', 'm35mms', 0.6215, T2), m35zams=('prod_c1_1124/m35zams/', 'm35zams', 0.6220, T2), t9l62=('prod_c4_1124/t9l62/', 't9l62', 0.674, T2), t9mlt=('prod_mlt_1125/t9mlt/', 't9l62', 0.674, T2))
SIG = 5.670374419e-5; MU_K = 1.66053906660e-24 / 1.380649e-16; TAB = {}
def table(f):
    g = [l for l in open(f) if re.match(r'#\s*\d', l)][0][1:].split(); nT, nD = int(g[0]), int(g[1]); lT0, dlT, lD0, dlD = map(float, g[2:6]); return nT, nD, lT0, dlT, lD0, dlD, np.loadtxt(f, comments='#').reshape(nT, nD)
def kap(t, rho, T):
    nT, nD, lT0, dlT, lD0, dlD, V = t; x = np.clip((np.log10(np.maximum(T, 1)) - lT0) / dlT, 0, nT - 1 - 1e-9); y = np.clip((np.log10(np.maximum(rho, 1e-30)) - lD0) / dlD, 0, nD - 1 - 1e-9); i = x.astype(int); j = y.astype(int); fx = x - i; fy = y - j
    return 10 ** ((1-fx)*(1-fy)*V[i, j] + fx*(1-fy)*V[i+1, j] + (1-fx)*fy*V[i, j+1] + fx*fy*V[i+1, j+1])
def asm(h, v):
    lg = np.asarray(h['mb_logical']); n3, n2, n1 = h['mb_data'][v][0].shape; P = np.empty((n3*(lg[:, 2].max()+1), n2*(lg[:, 1].max()+1), n1))
    for m in range(len(lg)): P[lg[m, 2]*n3:(lg[m, 2]+1)*n3, lg[m, 1]*n2:(lg[m, 1]+1)*n2] = h['mb_data'][v][m]
    return P
def one(a):
    s, n = a; fo = G + 'docfig_1125/p23/%s_%05d.npz' % (s, n)
    if os.path.exists(fo): return
    d, pre, mu, tf = R[s]
    if tf not in TAB: TAB[tf] = table(tf)
    h = bc.read_binary(G + d + 'bin/%s.hydro_w.%05d.bin' % (pre, n)); t = h['time'] / 86400; rho = asm(h, 'dens'); v = asm(h, 'velx'); T = (2/3.) * asm(h, 'eint') / rho * mu * MU_K; del h
    m = bc.read_binary(G + d + 'bin/%s.vet.%05d.bin' % (pre, n)); E = asm(m, 'vet_e'); F = asm(m, 'vet_f1'); del m
    gr = np.array([x.split()[:4] for x in open(G + d + pre + '.x1grid.txt') if not x.startswith('#') and len(x.split()) == 4], float); nr = rho.shape[2]; rc = gr[3:3+nr, 2]; dr = np.diff(gr[3:3+nr+1, 1])
    dt = kap(TAB[tf], rho, T) * rho * dr; del T
    tau = np.cumsum(dt[..., ::-1], -1)[..., ::-1] - 0.5 * dt; del dt
    above = tau < 2/3.; k = np.argmax(~above[..., ::-1], -1); ok = (~above).any(-1) & (k > 0); ib = np.clip(nr - 1 - k, 0, nr - 2)       # k = 0: the top cell itself is thick (photosphere off the top)
    g = lambda q, i: np.take_along_axis(q, i[..., None], -1)[..., 0]
    w = np.clip((g(tau, ib) - 2/3.) / np.maximum(g(tau, ib) - g(tau, ib+1), 1e-30), 0, 1); at = lambda q: g(q, ib) + (g(q, ib+1) - g(q, ib)) * w
    fc = F - 4/3. * E * v; f = g(fc, ib) + (g(fc, ib+1) - g(fc, ib)) * w
    np.savez_compressed(fo + '.tmp%d.npz' % os.getpid(), t=t, v=np.where(ok, at(v) / 1e5, np.nan).astype(np.float32), te=np.where(ok & (f > 0), (np.maximum(f, 1e-99) / SIG) ** 0.25, np.nan).astype(np.float32), r=np.where(ok, rc[ib] + (rc[ib+1] - rc[ib]) * w, np.nan).astype(np.float64), nnone=int((~ok).sum()), ntop=int((k == 0).sum()), ncol=ok.size); os.replace(fo + '.tmp%d.npz' % os.getpid(), fo); print(s, n, 't %.3f none %d of %d' % (t, (~ok).sum(), ok.size), flush=True)
if __name__ == '__main__':
    import multiprocessing as mp; s, r, w = sys.argv[1], sys.argv[2], int(sys.argv[3]); a, b = map(int, r.split('-'))
    with mp.get_context('spawn').Pool(w) as p: p.map(one, [(s, n) for n in range(a, b+1)], chunksize=1)
