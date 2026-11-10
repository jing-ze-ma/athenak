import sys, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires/code')
import csio
raw=csio.read_dump(sys.argv[1])
eos=csio.EOS('/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt')
mb=raw['mb_data']; rho=np.asarray(mb['dens'],float); T,p=eos.invert(rho,np.asarray(mb['eint'],float))
X,EXI,EET,dOm=csio.geometry(raw); rf,rc=csio.radii(raw)
print(rho.shape, X.shape, rf[0],rf[-1], raw['time'])
pc=np.median(p.reshape(-1,rho.shape[-1]),0)/1e6
i=np.argmin(abs(np.log10(pc)+2)); print('i at 1e-2 bar',i, rc[i])
Ti=T[...,i]
for ax in range(3): print(ax, (Ti*X[...,ax]*dOm).sum()/(Ti*dOm).sum())
print('Tmin max',Ti.min(),Ti.max())
for k in range(0,112,8): print(k, rc[k]/7.1492e9, pc[k], np.median(T[...,k]), T[...,k].min(), T[...,k].max())
