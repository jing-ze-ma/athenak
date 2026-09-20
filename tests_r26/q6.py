"""Q4: dt limiter census, floor atmosphere, grid-locked structure (settled dump 20)."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r26')
from three import (load, temp_of, pgas_of, weight, kappa, RS, TURN, ARAD, KB, MH,
                   VCEIL)  # noqa

CFL = 0.3
for run in ('lr1t', 'lr_f100t'):
    for idx in (20,):
        M, g = load(run, idx)
        rc, dr, wj = M['rc'], M['dr'], M['wj']
        n1, n2, n3 = M['n1'], M['n2'], M['n3']
        rho, vr, vt, vp, ei = (g['dens'], g['velx'], g['vely'], g['velz'], g['eint'])
        del g
        T = np.empty_like(rho)
        for a in range(0, n1, 16):
            b = min(a+16, n1)
            T[:, :, a:b] = temp_of(rho[:, :, a:b], ei[:, :, a:b])
        w = weight(np.log10(np.maximum(rho, 1e-300)))
        pg = pgas_of(rho, T)
        p = pg + w*ARAD*T**4/3.0
        beta = pg/np.maximum(p, 1e-300)
        gam1 = beta*(5./3.) + (1-beta)*(4./3.)
        cs = np.sqrt(gam1*p/np.maximum(rho, 1e-300))
        th = np.linspace(np.pi/4, 3*np.pi/4, n2+1)
        thc = 0.5*(th[1:]+th[:-1])
        dl2 = (rc[None, :]*(np.pi/2/n2))[None, :, :]*np.ones((1, n2, 1))
        dl3 = (rc[None, None, :]*np.sin(thc)[None, :, None]*(np.pi/2/n3))
        s1 = CFL*dr[None, None, :]/(np.abs(vr)+cs)
        s2 = CFL*dl2/(np.abs(vt)+cs)
        s3 = CFL*dl3/(np.abs(vp)+cs)
        sig = np.minimum(np.minimum(s1, s2), s3)
        print("\n" + "=" * 86)
        print("RUN %s dump %d t=%.3f turn : dt_hydro_min = %.4f s"
              % (run, idx, M['t']/TURN, sig.min()))
        k, j, i = np.unravel_index(sig.argmin(), sig.shape)
        which = ['radial', 'theta', 'phi'][int(np.argmin([s1[k, j, i], s2[k, j, i],
                                                          s3[k, j, i]]))]
        print("  limiting cell: i=%d r/R=%.4f (k=%d,j=%d) direction=%s  rho=%.3e T=%.3e "
              "cs=%.3e |v1|=%.3e |vt|=%.3e |vp|=%.3e  dr=%.3e beta=%.3f"
              % (i, rc[i]/RS, k, j, which, rho[k, j, i], T[k, j, i], cs[k, j, i],
                 abs(vr[k, j, i]), abs(vt[k, j, i]), abs(vp[k, j, i]), dr[i],
                 beta[k, j, i]))
        # what fraction of the 200 tightest cells sit where
        o = np.argsort(sig.ravel())[:200]
        kk, jj, ii = np.unravel_index(o, sig.shape)
        print("  the 200 tightest cells: r/R p10 %.3f med %.3f p90 %.3f ; "
              "cs-limited %.0f%% (|v|<cs) ; median cs %.2e |v| %.2e"
              % (np.percentile(rc[ii]/RS, 10), np.median(rc[ii]/RS),
                 np.percentile(rc[ii]/RS, 90),
                 100*np.mean(np.sqrt(vr**2+vt**2+vp**2).ravel()[o] < cs.ravel()[o]),
                 np.median(cs.ravel()[o]),
                 np.median(np.sqrt(vr**2+vt**2+vp**2).ravel()[o])))
        # per-radius minimum dt
        smin = sig.min(axis=(0, 1))
        kk = np.argsort(smin)[:8]
        print("  tightest radii: " + ", ".join("i=%d r/R=%.3f dt=%.3f"
                                               % (q, rc[q]/RS, smin[q]) for q in kk))
        # floor atmosphere
        top = rc > 1.05*RS
        print("  FLOOR ATMOSPHERE (r>1.05R): rho min %.3e med %.3e max %.3e ; "
              "T min %.3e med %.3e max %.3e ; frac with rho<1e-11 %.3f ; "
              "frac at T<6400 K %.4f ; frac with w>0 (LTE-tapered) %.3f"
              % (rho[:, :, top].min(), np.median(rho[:, :, top]), rho[:, :, top].max(),
                 T[:, :, top].min(), np.median(T[:, :, top]), T[:, :, top].max(),
                 (rho[:, :, top] < 1e-11).mean(), (T[:, :, top] < 6400).mean(),
                 (w[:, :, top] > 0).mean()))
        vmag = np.sqrt(vr**2+vt**2+vp**2)
        sel = vmag > 0.99*VCEIL
        if sel.any():
            print("  CEILING CELLS: %d ; r/R med %.3f p10 %.3f p90 %.3f ; rho med %.2e ;"
                  " T med %.2e ; beta med %.3f"
                  % (sel.sum(), np.median(rc[np.where(sel)[2]]/RS),
                     np.percentile(rc[np.where(sel)[2]]/RS, 10),
                     np.percentile(rc[np.where(sel)[2]]/RS, 90),
                     np.median(rho[sel]), np.median(T[sel]), np.median(beta[sel])))
        # grid-locked structure: Nyquist share on shells
        print("  GRID-LOCKED POWER (share of horizontal variance at 2-cell theta / phi /"
              " <=4 cells):")
        print("   %5s %7s | %-24s %-24s" % ("i", "r/R", "rho", "v_r"))
        for q in [2, 10, 21, 35, 52, 70, 85, 100, 120]:
            out = []
            for F in (rho[:, :, q], vr[:, :, q]):
                A = np.fft.fft2(F - F.mean())
                P = np.abs(A)**2
                P[0, 0] = 0.0
                tot = max(P.sum(), 1e-300)
                k3 = np.fft.fftfreq(n3)*n3
                k2 = np.fft.fftfreq(n2)*n2
                K3, K2 = np.meshgrid(k3, k2, indexing='ij')
                out.append((P[:, n2//2].sum()/tot, P[n3//2, :].sum()/tot,
                            P[(np.abs(K2) >= n2/4) | (np.abs(K3) >= n3/4)].sum()/tot))
            print("   %5d %7.4f | %7.4f %7.4f %7.4f    %7.4f %7.4f %7.4f"
                  % (q, rc[q]/RS, out[0][0], out[0][1], out[0][2],
                     out[1][0], out[1][1], out[1][2]))
