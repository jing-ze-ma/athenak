import sys, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc
fn=sys.argv[1]; phi0=float(sys.argv[2]); jj=int(sys.argv[3]) if len(sys.argv)>3 else 1
fd=bc.read_binary(fn); md=fd['mb_data']
for b in range(fd['n_mbs']):
    x3=fd['mb_x3v'][b]
    if x3[0]-1e-9 <= phi0 <= x3[-1]+1e-9:
        k=np.argmin(abs(x3-phi0)); r=fd['mb_x1v'][b]
        for i in list(range(50,75))+list(range(80,len(r),15)):
            print('%.4f rho %.4e vr %8.3f vp %8.3f T %.4g' % (r[i], md['dens'][b][k,jj,i], md['velx'][b][k,jj,i], md['velz'][b][k,jj,i], md['eint'][b][k,jj,i]/md['dens'][b][k,jj,i]*2/3))
        break
