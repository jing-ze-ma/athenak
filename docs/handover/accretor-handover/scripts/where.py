import sys, numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc
for fn in sys.argv[1:]:
    fd = bc.read_binary(fn)
    best = (-1, 0, 0, 0, 0, 0, 0, 0)
    for b in range(fd['n_mbs']):
        md = fd['mb_data']
        v = np.sqrt(md['velx'][b]**2 + md['velz'][b]**2)
        i = np.unravel_index(np.argmax(v), v.shape)
        if v[i] > best[0]:
            ee = md['eint'][b] if 'eint' in md else md[list(md.keys())[-1]][b]
            best = (v[i], fd['mb_x1v'][b][i[2]], fd['mb_x3v'][b][i[0]], md['dens'][b][i],
                    md['velx'][b][i], md['velz'][b][i], ee[i]/md['dens'][b][i]*(2/3), i[2])
    print(fn.split('/')[-1], 'max|v| %.4g at r %.4f phi %.3f rho %.4g vr %.4g vp %.4g P/rho %.4g i %d' % best)
