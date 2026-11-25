"""bsg concentric ripples near phi 30 deg just north of the equator: in which field, at which scale, in other dumps? read only"""
import sys, re, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d'); import bin_convert as bc
D='/resnick/groups/carnegie_poc/jingze/prod_c1_1124/bsg/'; RS=6.957e10
OP='/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt'
g=[l for l in open(OP) if re.match(r'#\s*\d',l)][0][1:].split(); nT,nD=int(g[0]),int(g[1]); lT0,dlT,lD0,dlD=map(float,g[2:6]); TAB=np.loadtxt(OP,comments='#').reshape(nT,nD)
def kap(rho,T):
    x=np.clip((np.log10(np.maximum(T,1))-lT0)/dlT,0,nT-1-1e-9); y=np.clip((np.log10(np.maximum(rho,1e-30))-lD0)/dlD,0,nD-1-1e-9); i=x.astype(int); j=y.astype(int); fx=x-i; fy=y-j
    return 10**((1-fx)*(1-fy)*TAB[i,j]+fx*(1-fy)*TAB[i+1,j]+(1-fx)*fy*TAB[i,j+1]+fx*fy*TAB[i+1,j+1])
gr=np.array([x.split()[:4] for x in open(D+'bsg3d.x1grid.txt') if not x.startswith('#') and len(x.split())==4],float); rc=gr[3:259,2]; rf=gr[3:260,1]; dr=np.diff(rf)
def asm(h,v):
    lg=np.asarray(h['mb_logical']); n3,n2,n1=h['mb_data'][v][0].shape; P=np.empty((n3*(lg[:,2].max()+1),n2*(lg[:,1].max()+1),n1))
    for m in range(len(lg)): P[lg[m,2]*n3:(lg[m,2]+1)*n3,lg[m,1]*n2:(lg[m,1]+1)*n2]=h['mb_data'][v][m]
    return P
def fields(n):
    h=bc.read_binary(D+'bin/bsg3d.hydro_w.%05d.bin'%n); t=h['time']/86400; rho=asm(h,'dens'); T=(2/3.)*asm(h,'eint')/rho*0.62*1.66053906660e-24/1.380649e-16; v=asm(h,'velx'); dt=str(np.asarray(h['mb_data']['dens'][0]).dtype); del h
    m=bc.read_binary(D+'bin/bsg3d.vet.%05d.bin'%n); E=asm(m,'vet_e'); F=asm(m,'vet_f1'); del m
    dt_=kap(rho,T)*rho*dr; tf=np.cumsum(dt_[...,::-1],-1)[...,::-1]; out={}
    tc=tf-0.5*dt_
    fc=F-4./3*E*v
    for lev in (2/3.,1.0):
        above=tf<lev; ib=np.argmax(above,-1)-1          # cell whose inner-face tau >= lev and next face < lev: crossing in cell ib
        ib=np.clip(ib,0,rc.size-2); a=np.take_along_axis
        ti=a(tf,ib[...,None],-1)[...,0]; to=a(tf,(ib+1)[...,None],-1)[...,0]; fr=np.clip((ti-lev)/np.maximum(ti-to,1e-30),0,1)
        rph=rf[ib]+fr*dr[ib]; out['r%.2f'%lev]=rph/RS
        # comoving flux interpolated linearly between cell centres to rph, and lab F
        def at(q):
            j=np.clip(np.searchsorted(rc,rph)-1,0,rc.size-2); w=(rph-rc[j])/(rc[j+1]-rc[j]); return a(q,j[...,None],-1)[...,0]*(1-w)+a(q,(j+1)[...,None],-1)[...,0]*w
        out['Fc%.2f'%lev]=at(fc); out['Fl%.2f'%lev]=at(F); out['E%.2f'%lev]=at(E); out['v%.2f'%lev]=at(v); out['rho%.2f'%lev]=at(np.log10(rho))
    # 'own tau=1 cell' reading used for the page and docs map: first cell whose face-integrated tau(from the top incl. the cell) < 1
    thin=tf<1.0; iown=np.where(thin.any(-1),np.argmax(thin,-1),rc.size-1); out['iown']=iown.astype(float); out['Fown']=np.take_along_axis(fc,iown[...,None],-1)[...,0]
    out['t']=t; out['Fshell']=F[...,int(np.argmin(abs(rc/RS-60)))]; out['Fs55']=F[...,int(np.argmin(abs(rc/RS-55)))]; out['Fs48']=fc[...,int(np.argmin(abs(rc/RS-48)))]; out['dtype']=dt; return out
if __name__=='__main__':
    for n in (int(x) for x in sys.argv[1:]):
        o=fields(n); np.savez_compressed('f_%05d.npz'%n,**{k:v for k,v in o.items() if k!='dtype'}); print(n,o['t'],o['dtype'],'rows done')
