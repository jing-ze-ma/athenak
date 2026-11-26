"""spec.py RUNDIR:NAME ... : lateral spectrum of the thin layer (cubed sphere 32/edge): per panel, 1-D FFT (Hann) of v2 along x2 and of v3 along x3
at the shells nearest 815 / 825 Rsun, power averaged over rows and the 6 panels; peak wavenumber k (per panel edge) and cells per wavelength N/k.
Dumps 60 / 70 / 75 d (rsg.hydro_w.000{12,14,15}.bin).  Same reader as conv_on.py."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.argv_save = sys.argv; sys.argv = ['x', 'cs', '.', 'x']
import conv_arm as C
sys.argv = sys.argv_save; RS = 6.957e10
for arg in sys.argv[1:]:
    D, nm = arg.split(':')
    for num in ('00012', '00014', '00015'):
        f = os.path.join(D, 'bin/rsg.hydro_w.%s.bin' % num)
        if not os.path.exists(f): continue
        t, cyc, bl = C.read(f); H, xf, xc, W, ang = C.cols(bl, ['dens', 'velx', 'vely', 'velz'])
        out = []
        for R in (815, 825):
            i = np.abs(xc - R*RS).argmin(); v2, v3, v1 = H['vely'][..., i], H['velz'][..., i], H['velx'][..., i]   # (6, N3, N2)
            N = v2.shape[-1]; w = np.hanning(N)
            def pw(a, ax):
                a = a - a.mean(axis=ax, keepdims=True); s = np.abs(np.fft.rfft(np.moveaxis(a, ax, -1)*w, axis=-1))**2; return s.reshape(-1, s.shape[-1]).mean(0)
            P = pw(v2, 2) + pw(v3, 1); Pr = pw(v1, 2) + pw(v1, 1); k = np.arange(P.size)
            kp = 1 + P[1:].argmax(); kr = 1 + Pr[1:].argmax(); hi = P[k >= N//4].sum()/P[1:].sum()
            out.append('%d Rsun: v_lat rms %.2f km/s, peak k %d (%.1f cells/wl), P(k>=%d)/P %.2f; v_r peak k %d (%.1f cells/wl) | Plat(k=1..10) %s' % (
                R, np.sqrt((v2**2 + v3**2).mean())/1e5, kp, N/kp, N//4, hi, kr, N/kr, ' '.join('%.2f' % x for x in P[1:11]/P[1:].max())))
        print('%-6s t %.1f d:\n   %s' % (nm, t/86400, '\n   '.join(out)), flush=True)
