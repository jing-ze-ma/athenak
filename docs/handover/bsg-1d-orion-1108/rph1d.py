import numpy as np,glob,sys
D=0.0669426; th=1e-10/6.15245e-5; nb=len(glob.glob('bsg.block*.out2.00000.tab'))
nmax=max(int(f.split('.')[-2]) for f in glob.glob('bsg.block0.out2.*.tab')); step=int(sys.argv[1]) if len(sys.argv)>1 else 20
for n in list(range(0,nmax+1,step))+([nmax] if nmax%step else []):
    a=np.vstack([np.loadtxt(f'bsg.block{b}.out2.{n:05d}.tab',comments='#',ndmin=2) for b in range(nb)])
    t=float(open(f'bsg.block0.out2.{n:05d}.tab').readline().split('time=')[1].split()[0])*D
    o=np.argsort(a[:,1]); f=a[o]; r=f[:,1]; sig=f[:,-3]+f[:,-2]
    rf=np.concatenate([[r[0]-(r[1]-r[0])/2],(r[1:]+r[:-1])/2,[r[-1]+(r[-1]-r[-2])/2]]); dr=np.diff(rf)
    tau=(sig*dr)[::-1].cumsum()[::-1]-.5*sig*dr
    i=np.nonzero(tau>1)[0][-1]
    rph=r[i]+(r[i+1]-r[i])*np.log(tau[i])/np.log(tau[i]/tau[i+1])
    j=np.where(f[:,2]<th)[0]
    print(f't={t:6.1f} d  r_photo(tau=1)={rph:6.2f}  r(rho=1e-10)={r[j[0]]:6.2f}')
