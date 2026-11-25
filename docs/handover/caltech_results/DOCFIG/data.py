import sys, re, numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d'); import bin_convert as bc
D='/resnick/groups/carnegie_poc/jingze/prod_c1_1124/bsg/'; RS=6.957e10
h=bc.read_binary(D+'bin/bsg3d.hydro_w.00114.bin'); lg=np.asarray(h['mb_logical']); print('t',h['time'],h['time']/86400)
def asm(v):
    n3,n2,n1=h['mb_data'][v][0].shape; P=np.empty((n3*(lg[:,2].max()+1),n2*(lg[:,1].max()+1),n1))
    for m in range(len(lg)): P[lg[m,2]*n3:(lg[m,2]+1)*n3,lg[m,1]*n2:(lg[m,1]+1)*n2]=h['mb_data'][v][m]
    return P
rho=asm('dens'); e=asm('eint'); v=asm('velx'); T=(2/3.)*e/rho*0.62*1.66053906660e-24/1.380649e-16
g=np.array([x.split()[:4] for x in open(D+'bsg3d.x1grid.txt') if not x.startswith('#') and len(x.split())==4],float); rc=g[3:259,2]; rf=g[3:260,1]
k=128; i=int(np.argmin(abs(rc/RS-40.64)))
c=np.load('/resnick/groups/carnegie_poc/jingze/pages_c1_1124/bsg/pipe/cache/view/int2_00114.npz')
np.savez('bsg114.npz',t=h['time'],rc=rc,rf=rf,T=T[k],rho=rho[k],vr=v[:,:,i],rs=rc[i],fr=c['fr_own'],rown=c['r_own'],th=np.linspace(np.pi/3,2*np.pi/3,257),ph=np.linspace(0,np.pi/3,257)); print('slice r',rc[i]/RS,'T range',T[k].min(),T[k].max(),'rho',rho[k].min(),rho[k].max(),'vr',v[:,:,i].min(),v[:,:,i].max(), 'fr',c['fr_own'].min(),c['fr_own'].max(), c['fr_own'].shape)
