import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc
fs=sorted(glob.glob(sys.argv[1]+'/bin/*.bin')); step=int(sys.argv[2]) if len(sys.argv)>2 else 1
for fn in fs[::step]+[fs[-1]]:
    fd=bc.read_binary(fn); md=fd['mb_data']; best=(1e30,)
    for b in range(fd['n_mbs']):
        r=fd['mb_x1v'][b]; x3=fd['mb_x3v'][b]; d=md['dens'][b]; e=md['eint'][b]
        c=np.sqrt(5/3*e*2/3/d)
        dr=np.gradient(r)[None,None,:]; dp=(x3[1]-x3[0])*r[None,None,:]
        t1=0.3*dr/(abs(md['velx'][b])+c); t3=0.3*dp/(abs(md['velz'][b])+c)
        dt=np.minimum(t1,t3); i=np.unravel_index(np.argmin(dt),dt.shape)
        if dt[i]<best[0]: best=(dt[i], 'r' if t1[i]<=t3[i] else 'phi', r[i[2]], x3[i[0]], d[i], md['velx'][b][i], md['velz'][b][i], c[i], e[i]/d[i]*2/3)
    print('t %.6f dt %.3g dir %s r %.4f phi %.3f rho %.3g vr %.4g vp %.4g c_ad %.4g P/rho %.4g'%((fd['time'],)+best))
