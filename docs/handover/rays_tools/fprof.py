import sys,types,numpy as np
sys.modules.setdefault('h5py', types.ModuleType('h5py')); sys.path.insert(0,'/orion/ptmp/jinma/clean_1115/int2_tests/ana'); import bin_convert as bc
exec(open('/orion/ptmp/jinma/clean_1115/stag2/zamsvphi/mir0.py').read().split('k=sys.argv[1]')[0])
run,k=sys.argv[1],sys.argv[2]; i0,i1=int(sys.argv[3]) if len(sys.argv)>3 else 320, int(sys.argv[4]) if len(sys.argv)>4 else 332
f3=full('%s/bin/m35zams.vet.%s.bin'%(run,k),'vet_f3')[:,:,i0:i1]; f1=full('%s/bin/m35zams.vet.%s.bin'%(run,k),'vet_f1')[:,:,i0:i1]
p=f3.sum(axis=(1,2))/np.abs(f1).sum(axis=(1,2)); K=len(p)
print(run,k,'i%d-%d F3/|F1|(k) every %d:'%(i0,i1-1,K//16),' '.join('%+.1e'%p[j] for j in range(0,K,K//16)),' | odd part max %.2e even(mean) %.2e'%(np.abs(p-p[::-1]).max()/2,p.mean()))
