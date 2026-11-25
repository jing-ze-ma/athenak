import numpy as np
from scipy.ndimage import uniform_filter
D='/resnick/groups/carnegie_poc/jingze/prod_c1_1124/bsg/'; RS=6.957e10
gr=np.array([x.split()[:4] for x in open(D+'bsg3d.x1grid.txt') if not x.startswith('#') and len(x.split())==4],float); rf=gr[3:260,1]/RS
ph=(np.arange(256)+.5)*60/256; lat=30-(np.arange(256)+.5)*60/256
def hp(q,n=9): return q-uniform_filter(q,n,mode='nearest')
for n in (114,113,110,104):
    z=np.load('f_%05d.npz'%n)
    for nm,(a,b,c,d) in (('phi 26-36 lat 4-11',(26,36,4,11)),('whole map',(0,60,-30,30))):
        s=np.ix_(np.where((ph>a)&(ph<b))[0],np.where((lat>c)&(lat<d))[0]); out=[]
        for lev in ('0.67','1.00'):
            x=np.interp(z['r'+lev],rf,np.arange(rf.size)); f=x%1.0     # position of the crossing inside its radial cell (0..1)
            for k in ('Fshell','Fs55','Fc'+lev,'Fown','Fs48'):
                q=hp(z[k].astype(float))/uniform_filter(z[k].astype(float),9,mode='nearest'); q=q[s].ravel(); ff=f[s].ravel()
                M=np.stack([np.cos(2*np.pi*ff),np.sin(2*np.pi*ff),np.ones(ff.size)],1); co,res,_,_=np.linalg.lstsq(M,q,rcond=None); r2=1-((q-M@co)**2).sum()/((q-q.mean())**2).sum()
                # null: same fit with the position field shifted by 40 cells (destroys the link)
                f2=np.roll(f,40,axis=0)[s].ravel(); M2=np.stack([np.cos(2*np.pi*f2),np.sin(2*np.pi*f2),np.ones(f2.size)],1); c2=np.linalg.lstsq(M2,q,rcond=None)[0]; n2=1-((q-M2@c2)**2).sum()/((q-q.mean())**2).sum()
                out.append('%s(tau %s) R2 %.2f (null %.2f, rms %.3f)'%(k,lev,r2,n2,q.std()))
        print(n,nm,' | '.join(out))
# number of cells crossed in the region
z=np.load('f_00114.npz'); s=np.ix_(np.where((ph>26)&(ph<36))[0],np.where((lat>4)&(lat<11))[0]); print('r(tau=2/3) in region: %.2f-%.2f Rsun = %.1f radial cells of %.3f Rsun'%(z['r0.67'][s].min(),z['r0.67'][s].max(),(z['r0.67'][s].max()-z['r0.67'][s].min())/0.15625,0.15625))
