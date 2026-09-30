"""emis_prt.py: angle-dependent emission of every column with the pRT R1000 premixed table.

  python emis_prt.py <arm w1x|w10x> <table 1x_eq|10x_eq|1x_nodiss> [nproc=16]
Streams columns over a process pool (fork: the table is shared read-only).  Per column: k from
the table (bilinear in log T, log p of log10 k; clamped), dtau = k rho dr (the run's rho and
radial grid), source r^2 B_lambda(T) at bin centres, formal solution (rtcore.formal32, float32) at the
8 Gauss nodes and at the column's mu toward the observer for each of 72 phases.
Output out/prt_<arm>_<table>.npz:
  Pang (72, nfreq)   disk flux x d^2 [erg/s/cm] per unit wavelength, angle-dependent
  Fcol (ncol, nfreq) float32: top-face upward flux spectrum (8-pt Gauss) [erg/s/cm2/cm]
  lam, lam_edges [um], phi, dOm, rtop, X
"""
import sys
import numpy as np
from multiprocessing import Pool
import rtcore as rc

G = {}


def init(arm, tab):
    st = np.load(rc.SD + 'state/%s.npz' % arm)
    t = np.load(rc.SD + 'tables/prt_%s.npz' % tab)
    G.update(T=st['T'], p=st['p']/1e6, rho=st['rho'], rf=st['rf'], X=st['X'], dOm=st['dOm'],
             lT=t['lT'], lP=t['lP'], lk=t['lk'], gw=t['gw'].astype(np.float64),
             lam=t['lam'], le=t['lam_edges'])
    phi, obs = rc.phases(72)
    G['phi'] = phi
    G['mup'] = G['X'] @ obs.T


def interp(x, grid):
    f = (np.clip(x, grid[0], grid[-1]) - grid[0])/(grid[1] - grid[0])
    i = np.clip(f.astype(int), 0, len(grid) - 2)
    return i, (f - i).astype(np.float32)


def column(q):
    lk = G['lk']
    iT, fT = interp(np.log10(G['T'][q]), G['lT'])
    iP, fP = interp(np.log10(G['p'][q]), G['lP'])
    a, b = fT[:, None, None], fP[:, None, None]
    L = ((1 - a)*((1 - b)*lk[iT, iP] + b*lk[iT, iP + 1])
         + a*((1 - b)*lk[iT + 1, iP] + b*lk[iT + 1, iP + 1]))       # (nlay, nf, ng)
    rf = G['rf']
    dz = np.diff(rf)
    rcen = 0.5*(rf[1:] + rf[:-1])
    dtau = (10.0**L)*(G['rho'][q]*dz)[:, None, None].astype(np.float32)
    lamc = G['lam']*1e-4
    S = (rc.planck_lam(lamc[None, :], G['T'][q][:, None])*rcen[:, None]**2)[:, :, None]
    gx8, gw8 = rc.gauss01(8)
    vis = G['mup'][q] > 0
    mus = np.concatenate([gx8, G['mup'][q, vis]])
    J = rc.formal32(dtau, rc.face_source(S, dtau), mus)                  # (nf, ng, nmu)
    Jg = np.einsum('fgm,g->fm', J, G['gw'])
    F = 2*np.pi*(Jg[:, :8]*(gx8*gw8)).sum(1)/rf[-1]**2
    return q, F.astype(np.float32), vis, Jg[:, 8:]*(G['mup'][q, vis]*G['dOm'][q])


def chunk(qs):
    nf = len(G['lam'])
    P = np.zeros((72, nf))
    out = []
    for q in qs:
        q, F, vis, c = column(q)
        P[vis] += c.T
        out.append((q, F))
    return P, out


if __name__ == '__main__':
    arm, tab = sys.argv[1], sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    init(arm, tab)
    ncol = G['T'].shape[0]
    nf = len(G['lam'])
    Pang = np.zeros((72, nf))
    Fcol = np.zeros((ncol, nf), np.float32)
    chunks = [list(range(s, min(s + 32, ncol))) for s in range(0, ncol, 32)]
    # Robust to lost workers (09-30: 5 workers died with SIGBUS when the shared 70 GB
    # claude-work.slice hit MemoryMax; a dead Pool worker silently loses its task, so every
    # chunk is waited for with a timeout and resubmitted, and completeness is asserted).
    done = np.zeros(ncol, bool)
    with Pool(nproc) as pool:
        pend = {i: pool.apply_async(chunk, (c,)) for i, c in enumerate(chunks)}
        tries = {i: 1 for i in pend}
        n = 0
        while pend:
            for i in list(pend):
                try:
                    P, out = pend[i].get(timeout=1200)
                except Exception as ex:
                    if tries[i] >= 4:
                        raise
                    print('chunk %d lost (%s), resubmitting' % (i, type(ex).__name__), flush=True)
                    pend[i] = pool.apply_async(chunk, (chunks[i],))
                    tries[i] += 1
                    continue
                del pend[i]
                Pang += P
                for q, F in out:
                    Fcol[q] = F
                    done[q] = True
                n += 1
                if n % 20 == 0:
                    print('chunk %d / %d' % (n, len(chunks)), flush=True)
    assert done.all(), 'columns missing: %d' % (~done).sum()
    print('all %d columns done; resubmitted chunks: %d' % (ncol, sum(t > 1 for t in tries.values())))
    rtop = G['rf'][-1]
    dl = np.diff(G['le'])*1e-4 if G['le'][1] > G['le'][0] else -np.diff(G['le'])*1e-4
    np.savez(rc.SD + 'out/prt_%s_%s.npz' % (arm, tab), Pang=Pang, Fcol=Fcol, lam=G['lam'],
             lam_edges=G['le'], phi=G['phi'], dOm=G['dOm'], rtop=rtop, X=G['X'])
    L = (Fcol.astype(float)*np.abs(dl)).sum(1) @ G['dOm']*rtop**2
    print(arm, tab, 'L_IR(0.3-28 um) = %.4e erg/s' % L)
