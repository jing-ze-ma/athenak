import sys, os, numpy as np
import bin_convert as bc
B='/raven/ptmp/jinma/bsgnv_1107/long3/bin/bsg3d.'
Rs=6.955e10
for n in sys.argv[1:]:
    n=int(n)
    m=bc.read_binary(B+'m1.%05d.bin'%n); v=bc.read_binary(B+'m1_vet.%05d.bin'%n)
    h=bc.read_binary(B+'hydro_w.%05d.bin'%n)
    fn=B+'m1_face.%05d.bin'%n; fc=bc.read_binary(fn) if os.path.exists(fn) else None
    if n==int(sys.argv[1]) and fc is not None: print('face vars',fc['var_names'])
    x1=np.array(m['mb_x1v'][0]); 
    def sh(d,k): return np.asarray(d['mb_data'][k]).mean(axis=(0,1,2))
    E=sh(m,'m1_e'); F=sh(m,'m1_f1'); frr=sh(v,'m1_d11'); rho=sh(h,'dens')
    Ff = sh(fc,'m1_fx1') if fc is not None and 'm1_fx1' in fc['var_names'] else None
    print('dump %d  t = %.2f d'%(n,m['time']/86400))
    for i in range(0,256,8):
        r=x1[i]/Rs
        if r<48: continue
        s='  r %6.2f rho %8.2e F/cE(cell,capped) %8.2e f_rr %6.3f'%(r,rho[i],F[i]/E[i],frr[i])
        if Ff is not None: s+='  Fface/cE %8.2e'%(Ff[i]/E[i])
        print(s)
