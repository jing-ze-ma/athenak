import sys, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc
for fn in sys.argv[1:]:
    fd = bc.read_binary(fn)
    md = fd['mb_data']
    best=(-1,)
    for b in range(fd['n_mbs']):
        r = fd['mb_x1v'][b]
        v = np.sqrt(md['velx'][b]**2+md['velz'][b]**2)
        m = r[None,None,:] < 4.06
        vv = np.where(m, v, -1)
        i = np.unravel_index(np.argmax(vv), vv.shape)
        if vv[i] > best[0]:
            best=(vv[i], r[i[2]], fd['mb_x3v'][b][i[0]], md['dens'][b][i], md['velx'][b][i], md['eint'][b][i]/md['dens'][b][i]*2/3, b, i)
    print(fn.split('.')[-2], 't %.5f'%fd['time'], 'env max|v| %.4g r %.4f phi %.3f rho %.4g vr %.4g P/rho %.4g'%best[:6], best[6:])
