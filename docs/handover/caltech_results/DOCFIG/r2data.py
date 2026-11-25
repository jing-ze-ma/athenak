import sys, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d'); import bin_convert as bc
D='/resnick/groups/carnegie_poc/jingze/prod_c1_1124/bsg/'; RS=6.957e10
h=bc.read_binary(D+'bin/bsg3d.hydro_w.00114.bin'); lg=np.asarray(h['mb_logical'])
def asm(v):
    n3,n2,n1=h['mb_data'][v][0].shape; P=np.empty((n3*(lg[:,2].max()+1),n2*(lg[:,1].max()+1),n1))
    for m in range(len(lg)): P[lg[m,2]*n3:(lg[m,2]+1)*n3,lg[m,1]*n2:(lg[m,1]+1)*n2]=h['mb_data'][v][m]
    return P
rho=asm('dens'); e=asm('eint'); v=asm('velx'); T=(2/3.)*e/rho*0.62*1.66053906660e-24/1.380649e-16
g=np.array([x.split()[:4] for x in open(D+'bsg3d.x1grid.txt') if not x.startswith('#') and len(x.split())==4],float); rc=g[3:259,2]; rf=g[3:260,1]
th=np.linspace(np.pi/3,2*np.pi/3,257); w=np.cos(th[:-1])-np.cos(th[1:]); w=w/w.sum()/256
sm=lambda q: np.einsum('kji,j->i',q,w)
rm,Tm=sm(rho),sm(T); vr2=np.sqrt(sm((v-sm(v))**2)); dr=np.sqrt(sm((rho/rm-1)**2)); dT=np.sqrt(sm((T/Tm-1)**2)); r=rc/RS
fl=6.15245e-14
for i in range(0,256,6): print('r %.2f rms v %.2f km/s drho/rho rms %.3f dT/T %.4f floor frac %.2f rho %.2e'%(r[i],vr2[i]/1e5,dr[i],dT[i],(rho[:,:,i]<=1.05*fl).mean(),rm[i]))
np.savez('r2_shell.npz',r=r,vr2=vr2,dr=dr,dT=dT)
np.savez_compressed('r2_vol.npz',rho=rho[:,:,:200].astype(np.float32),T=T[:,:,:200].astype(np.float32),v=v[:,:,:200].astype(np.float32),rm=rm,Tm=Tm,r=r,rf=rf/RS,fl=fl,t=h['time'])
