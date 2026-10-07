import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc
# dt-limiting cell (cfl 0.3, gamma 5/3), and mass outside R_acc (ambient + stream) per dump
for fn in sorted(glob.glob(sys.argv[1]+'/bin/*.bin')):
    fd=bc.read_binary(fn); md=fd['mb_data']; best=(1e30,); mo=0.0; mfl=0
    for b in range(fd['n_mbs']):
        r=fd['mb_x1v'][b]; x3=fd['mb_x3v'][b]; d=md['dens'][b]
        c=np.sqrt(5/3*md['eint'][b]*2/3/d)
        dr=np.gradient(r)[None,None,:]; dp=(x3[1]-x3[0])*r[None,None,:]
        dt=0.3*np.minimum(dr/(abs(md['velx'][b])+c), dp/(abs(md['velz'][b])+c))
        i=np.unravel_index(np.argmin(dt),dt.shape)
        if dt[i]<best[0]: best=(dt[i], r[i[2]], x3[i[0]], d[i], md['velx'][b][i], md['velz'][b][i], c[i])
        out=r[None,None,:]>4.06
        vol=(r**2*dr[0,0])[None,None,:]*np.ones_like(d)
        mo+=np.sum(np.where(out,d*vol,0)); mfl+=np.sum(out&(d<=1.0001e-7))
    print('t %.4f dt_min %.3g at r %.4f phi %.3f rho %.3g vr %.4g vp %.4g c %.4g | M(r>R_acc) %.4g (r^2 dr dphi dtheta omitted), floor cells %d'%((fd['time'],)+best+(mo,mfl)))
