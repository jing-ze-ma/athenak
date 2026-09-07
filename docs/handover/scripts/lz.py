import sys, glob, numpy as np
sys.path.insert(0,'/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'); import dhjcs
C=[-0.068392,-2.191487,2.464818,-1.366698]
def stretch_r(r,r0,r1):
    xi=(r-r0)/(r1-r0); u=xi.copy(); xik=xi.copy()
    for k in range(1,5): u=u+C[k-1]*xik*(1-xi); xik=xik*xi
    return r0+(r1-r0)*u
def lz(path,kind):
    raw=dhjcs.bin_convert.read_binary(path); mb=raw['mb_data']; g=np.asarray(raw['mb_geometry']); rho=np.asarray(mb['dens'])
    nmb,nk,nj,ni=rho.shape; r0,r1=raw['x1min'],raw['x1max']
    rf=stretch_r(g[0,0]+(g[0,1]-g[0,0])*np.arange(ni+1)/ni, r0,r1); rc=0.5*(rf[1:]+rf[:-1]); dV_r=(rf[1:]**3-rf[:-1]**3)/3
    if kind=='cs':
        geo=dhjcs.cs_geometry(raw); u,v=dhjcs.winds(raw,geo); lat=geo['lat']
        dxi=np.pi/4*(g[0,3]-g[0,2])/nj; deta=np.pi/4*(g[0,5]-g[0,4])/nk
        A=np.zeros((nmb,nk,nj))
        for m in range(nmb):
            xi=np.pi/4*(g[m,2]+(g[m,3]-g[m,2])*(np.arange(nj)+0.5)/nj); eta=np.pi/4*(g[m,4]+(g[m,5]-g[m,4])*(np.arange(nk)+0.5)/nk)
            XI,ETA=np.meshgrid(xi,eta); X=np.tan(XI); Y=np.tan(ETA); d=np.sqrt(1+X*X+Y*Y); A[m]=(1+X*X)*(1+Y*Y)/d**3*dxi*deta
    else:
        geo=dhjcs.sp_geometry(raw); u,v=dhjcs.sp_winds(raw); lat=geo['lat']
        A=np.zeros((nmb,nk,nj))
        for m in range(nmb):
            thf=g[m,2]+(g[m,3]-g[m,2])*np.arange(nj+1)/nj; xi=thf/np.pi; thf=np.pi/2*(1+np.sinh(3.0*(2*xi-1))/np.sinh(3.0))
            dph=(g[m,5]-g[m,4])/nk; A[m]=np.tile((np.cos(thf[:-1])-np.cos(thf[1:]))*dph,(nk,1))
    dV=A[...,None]*dV_r[None,None,None,:]
    M=(rho*dV).sum(); L=(rho*u*rc[None,None,None,:]*np.cos(lat)[...,None]*dV).sum(); Lr=(rho*(2.06e-5)*rc[None,None,None,:]**2*np.cos(lat)[...,None]**2*dV).sum()
    return raw['time']/3.05e5, M, L, Lr, 4/3*np.pi*(r1**3-r0**3), dV.sum()
B='/viper/u2/jinma/ATHENAK/bench/'
for name,kind,fs in [('cs hydro','cs',[B+'cs_prod_hyd_rot/bin/dhj.hydro_w.%05d.bin'%i for i in [0,1,5,10,25,50,100,135]]),
                     ('cs MHD bcfix','cs',[B+'cs_inflow_bcfix/bin/dhj.mhd_w_bcc.%05d.bin'%i for i in [8,14,20]]),
                     ('cs MHD prod (new)','cs',sorted(glob.glob(B+'cs_mhd_prod/bin/*.bin'))),
                     ('sp hydro','sp',[B+'sp_dhj_hyd/bin/dhj.hydro_w.%05d.bin'%i for i in [0,1,5,10,25,50,83]])]:
    print('==',name); L0=None
    for f in fs:
        try: t,M,L,Lr,Vex,V=lz(f,kind)
        except Exception as e: print('  skip',f,e); continue
        if L0 is None: L0=L; M0=M; print(f'   volume check: sum dV / exact = {V/Vex:.6f}')
        print(f'   rot {t:6.1f}  M/M0-1 = {M/M0-1:+.2e}   Lz(rot frame) = {L:+.4e}  Lz/(Omega-frame L) = {L/Lr:+.4e}   Lz - Lz0 relative to Omega-frame L = {(L-L0)/Lr:+.3e}')
