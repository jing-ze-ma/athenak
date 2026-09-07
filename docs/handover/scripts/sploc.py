import sys, glob, numpy as np
sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/run'); import bin_convert
for f in sys.argv[1:]:
    raw=bin_convert.read_binary(f); mb=raw['mb_data']; t=raw['time']
    if 'first' not in globals():
        first=1; print('keys', [k for k in raw.keys() if k!='mb_data'])
    rho=np.asarray(mb['dens']); vr=np.asarray(mb['velx']); vt=np.asarray(mb['vely']); vp=np.asarray(mb['velz'])
    b=np.sqrt(sum(np.asarray(mb[k])**2 for k in ['bcc1','bcc2','bcc3']))
    geo=None
    for key in ['mb_geometry','mb_logical']:
        if key in raw: geo=np.asarray(raw[key]); print(key, geo.shape) if first==1 else None
    first=2
    b0=b[:,:,:,0]; idx=np.unravel_index(np.nanargmax(b0), b0.shape)
    print(f'{f}: t={t:.3e} rot={t/3.05e5:.2f}')
    print(f'  Bmax(i=0)={b0.max():.3e} at (m,k,j)={idx}; Bmax(all)={np.nanmax(b):.3e} at {np.unravel_index(np.nanargmax(b), b.shape)}')
    print(f'  rho min={np.nanmin(rho):.3e} at {np.unravel_index(np.nanargmin(rho), rho.shape)}; |v| max={np.nanmax(np.sqrt(vr**2+vt**2+vp**2))/1e5:.2f} km/s at {np.unravel_index(np.nanargmax(vr**2+vt**2+vp**2), vr.shape)}; vr min={np.nanmin(vr)/1e5:.3f} at {np.unravel_index(np.nanargmin(vr), vr.shape)}')
    if 'mb_geometry' in raw:
        g=np.asarray(raw['mb_geometry']); m=idx[0]; print(f'  MB {m} geometry x1[{g[m,0]:.3e},{g[m,1]:.3e}] x2[{g[m,2]:.4f},{g[m,3]:.4f}] x3[{g[m,4]:.4f},{g[m,5]:.4f}]')
    # radial profile of max B over all cells at each i, first 8 i and mean
    print('  max|B| by i (0..7):', ' '.join(f'{np.nanmax(b[:,:,:,i]):.2e}' for i in range(8)))
    print('  count |B|>1 kG by i (0..7):', ' '.join(f'{int((b[:,:,:,i]>1e3).sum())}' for i in range(8)))
