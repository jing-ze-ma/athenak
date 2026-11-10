import sys, numpy as np, csio
eos = csio.EOS('/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt')
for f in sys.argv[1:]:
    r = csio.read_dump(f); mb = r['mb_data']; rho = np.asarray(mb['dens'], float)
    T, p = eos.invert(rho, np.asarray(mb['eint'], float)); pb = p/1e6
    X, EXI, EET, _ = csio.geometry(r)
    vr, v2, v3 = [np.asarray(mb[k], float) for k in ('velx', 'vely', 'velz')]
    m = (pb > 1e-5) & (pb < 1e-3)
    # zonal velocity (about z): e_phi = z x rhat / |..|
    ez = np.array([0, 0, 1.0]); ph = np.cross(ez, X); ph /= np.maximum(np.linalg.norm(ph, axis=-1), 1e-12)[..., None]
    vphi = v2*np.sum(EXI*ph, -1)[..., None] + v3*np.sum(EET*ph, -1)[..., None]
    vh = np.sqrt(v2**2 + v3**2)
    print(f.split('/')[-3], f[-9:-4], 't/P %.2f' % (r['time']/235356.5), 'p 1e-5..1e-3 bar: <vphi> %.2f km/s  rms vh %.2f  max vh %.2f  rms vr %.3f  max|vr| %.2f  max T %.0f' % (
        vphi[m].mean()/1e5, np.sqrt((vh[m]**2).mean())/1e5, vh[m].max()/1e5, np.sqrt((vr[m]**2).mean())/1e5, np.abs(vr[m]).max()/1e5, T[m].max()))
