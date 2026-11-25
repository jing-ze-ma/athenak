import glob, numpy as np
RP=9.554e11
def load(run):
    S=[dict(np.load(f)) for f in sorted(glob.glob('m35c/%s/s_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); return S
A=load('m35c'); B=load('c1'); r=A[0]['r']/RP; ir=lambda x:int(np.argmin(abs(r-x))); ta=np.array([float(z['t']) for z in A]); tb=np.array([float(z['t']) for z in B])
bands={'A 9-14 (same cells: 128 grid 9-14 per wavelength / 256 grid k 11-12 -> 11 cells)':(9,14),'B 19-27 (l about 400)':(19,27),'lo 2-6':(2,6),'C 30-45':(30,45)}
print('dumps m35c',len(A),ta[0],ta[-1],'c1',len(B),tb[0],tb[-1],'count share of axis modes in band:',{k:round(float(A[0]['cax'][a:b+1].sum()/A[0]['cnt'][a:b+1].sum()),3) for k,(a,b) in bands.items()},{k:round(float(B[0]['cax'][a:b+1].sum()/B[0]['cnt'][a:b+1].sum()),3) for k,(a,b) in bands.items()})
def rms(z,i,a,b): return np.sqrt(z['sp'][a:b+1,i].sum())/1e5
def ax(z,i,a,b): return z['sa'][a:b+1,i].sum()/max(z['sp'][a:b+1,i].sum(),1e-99)
def near(S,T,t): return S[int(np.argmin(abs(T-t)))]
for x in (0.72,0.74,0.76):
    i=ir(x); print('=== r = %.2f R_ph,IC'%x)
    for tt in (1.0e5,1.4e5,1.8e5,2.16e5,2.6e5,3.0e5,3.46e5):
        za,zb=near(A,ta,tt),near(B,tb,tt); out='t %.2e (c1 %.2e, m35c %.2e): '%(tt,float(zb['t']),float(za['t']))
        for nm,(a,b) in bands.items():
            out+='[%s] c1 %.2e m35c %.2e km/s axis %.2f/%.2f; '%(nm.split()[0],rms(zb,i,a,b),rms(za,i,a,b),ax(zb,i,a,b),ax(za,i,a,b))
        pk=lambda z:int(6+np.argmax(z['sp'][6:46,i])); out+='peak k (6-45): c1 %d m35c %d | total rms c1 %.3e m35c %.3e'%(pk(zb),pk(za),zb['rms'][i]/1e5,za['rms'][i]/1e5); print('  '+out)
    for nm,(a,b) in bands.items():
        for S,T,lab in ((B,tb,'c1'),(A,ta,'m35c')):
            m=(T>=1.0e5)&(T<=3.5e5); y=np.array([rms(z,i,a,b) for z in S])[m]; tt=T[m]; p=np.polyfit(tt,np.log(y+1e-30),1); o=[]
            for u,v in ((1.0e5,2.0e5),(2.0e5,2.6e5),(2.6e5,3.5e5)):
                q=(tt>=u)&(tt<=v)
                if q.sum()>3: pp=np.polyfit(tt[q],np.log(y[q]+1e-30),1); o.append('%.2e-%.2e: %.2e s'%(u,v,1/pp[0]) if pp[0]!=0 else '')
            print('   growth e-fold band %s %s: %s'%(nm.split()[0],lab,' | '.join(o)))
# full spectrum at the last common time, at 0.74
i=ir(0.74); za,zb=near(A,ta,3.46e5),near(B,tb,3.46e5); print('spectrum at 0.74, t 3.46e5 (share of rms^2 per k bin, percent): c1',[(k,round(100*float(zb['sp'][k,i]/zb['sp'][1:,i].sum()),1)) for k in range(1,46) if zb['sp'][k,i]/zb['sp'][1:,i].sum()>0.02]); print('  m35c',[(k,round(100*float(za['sp'][k,i]/za['sp'][1:,i].sum()),1)) for k in range(1,46) if za['sp'][k,i]/za['sp'][1:,i].sum()>0.02])
