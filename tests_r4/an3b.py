import os, sys, numpy as np
kB=1.380649e-16; mH=1.6726219e-24
cases={}
for row in open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','tests_r2','thick','cases.txt')):
    if row.startswith('#') or not row.strip(): continue
    f=row.split(); cases[f[0]]=(float(f[4]),float(f[5]))   # x1max, kappa
print("%-16s %5s %9s %9s %9s %9s"%("case","N","min","max","median","max|1-x|"))
for tag,(x1max,kap) in cases.items():
    d=np.loadtxt(sys.argv[1]+"/"+tag+"/mltfaces.txt")
    r,T,p=d[:,1],d[:,2],d[:,3]; freq,f2s=d[:,10],d[:,17]
    dr=(x1max-1.0e12)/128.0
    rho=p*mH/(kB*T)            # mu = 1, ideal
    dtau=kap*rho*dr
    ok=(dtau>0.3)&(freq!=0)
    ok[:6]=False; ok[-6:]=False   # the wall face and the top ghost mirror
    rr=(f2s/freq)[ok]
    print("%-16s %5d %9.5f %9.5f %9.5f %9.5f"%(tag,ok.sum(),rr.min(),rr.max(),
                                               np.median(rr),np.abs(rr-1).max()))
