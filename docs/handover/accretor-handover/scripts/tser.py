import sys, glob, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc
# per dump: max T and max v over r in [4.07, 4.6] (all phi), and where; T at r index for phi of max
for fn in sorted(glob.glob(sys.argv[1]+'/bin/*.bin')):
    fd=bc.read_binary(fn); md=fd['mb_data']; best=(-1,)
    for b in range(fd['n_mbs']):
        r=fd['mb_x1v'][b]; T=md['eint'][b]/md['dens'][b]*2/3
        msk=(r[None,None,:]>4.0)&(r[None,None,:]<4.6)
        TT=np.where(msk,T,-1); i=np.unravel_index(np.argmax(TT),TT.shape)
        if TT[i]>best[0]: best=(TT[i], r[i[2]], fd['mb_x3v'][b][i[0]], md['dens'][b][i], md['velx'][b][i], i[1])
    print('t %.4f maxT %.4g at r %.4f phi %.3f rho %.3g vr %.4g j %d'%((fd['time'],)+best))
