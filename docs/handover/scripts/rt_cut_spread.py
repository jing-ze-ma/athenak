import sys, numpy as np
S='/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'
sys.path.insert(0,S); sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/vis/python')
import dhjcs, lz
eos=dhjcs.EOS(S+'/eostab/eos_table.txt')
def freedman(T,p,met=0.0):
    T1=np.clip(T,75,4000); p1=np.clip(p,1,3e8); lT=np.log10(T1); lp=np.log10(p1)
    c1,c2,c3,c4,c5,c6,c7=10.602,2.882,6.09e-15,2.954,-2.526,0.843,-5.490
    lo=T1<800
    c8=np.where(lo,-14.051,82.241);c9=np.where(lo,3.055,-55.456);c10=np.where(lo,0.024,8.754)
    c11=np.where(lo,1.877,0.7048);c12=np.where(lo,-0.445,-0.0414);c13=0.8321
    lkl=c1*np.arctan(lT-c2)-c3/(lp+c4)*np.exp((lT-c5)**2)+c6*met+c7
    lkh=c8+c9*lT+c10*lT**2+lp*(c11+c12*lT)+c13*met*(0.5+np.arctan((lT-2.5)/0.2)/np.pi)
    return 10**lkl+10**lkh
path=sys.argv[1]; kind=sys.argv[2]
d=dhjcs.load(path,eos,kind); raw=d['raw']; g=np.asarray(raw['mb_geometry']); rho=d['rho']; T=d['T']; p=d['p']
nmb,nk,nj,ni=rho.shape; r0,r1=raw['x1min'],raw['x1max']
rf=lz.stretch_r(g[0,0]+(g[0,1]-g[0,0])*np.arange(ni+1)/ni, r0,r1); dr=np.diff(rf)
kr=freedman(T,p)   # p in dyn/cm2? check
dtau=kr*rho*dr[None,None,None,:]
tau=np.cumsum(dtau[...,::-1],axis=-1)[...,::-1]   # tau at the BOTTOM face of each cell... approx: tau(top face of cell i)=sum over cells above
tau_top=tau-dtau   # tau at the top face of cell i
# RT bottom face = first face (from the bottom) with w<1, i.e. tau_face < 300; face i is the bottom face of cell i => tau_face(i)=tau(i)
icut=np.argmax(tau<300.0,axis=-1)   # first cell index (bottom up) whose bottom face has tau<300
old=np.argmax(p<10.0e6,axis=-1)
print('rot %.1f  %s  ni=%d'%(d['rot'],kind,ni))
for name,ic in [('tau_R=300 cut',icut),('10 bar cut',old)]:
    L=ni-ic
    print('%-14s bottom cell: min %d  median %d  max %d | chain length: min %d mean %.1f max %d | warp-loss 1-mean/max = %.2f'%(name,ic.min(),np.median(ic),ic.max(),L.min(),L.mean(),L.max(),1-L.mean()/L.max()))
    pc=np.take_along_axis(p,ic[...,None],axis=-1)[...,0]/1e6
    print('   pressure at that face: min %.3g  median %.3g  max %.3g bar'%(pc.min(),np.median(pc),pc.max()))
ov=(np.argmax(tau<30.0,axis=-1)-icut); print('overlap width cells: min %d mean %.1f max %d'%(ov.min(),ov.mean(),ov.max()))
