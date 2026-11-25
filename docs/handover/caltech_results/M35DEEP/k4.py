import glob, numpy as np
RP=9.554e11
def load(run):
    S=[dict(np.load(f)) for f in sorted(glob.glob('m35c/%s/s_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); return S
A=load('m35c'); B=load('c1'); r=A[0]['r']/RP; ir=lambda x:int(np.argmin(abs(r-x)))
for x in (0.72,0.74,0.76):
    i=ir(x); print('r',x)
    for tt in (1.0e5,1.4e5,1.8e5,2.16e5,2.6e5,3.0e5,3.46e5):
        o='  t %.2e '%tt
        for lab,S in (('c1',B),('m35c',A)):
            z=min(S,key=lambda q: abs(float(q['t'])-tt)); sp=z['sp'][:,i]; k=np.arange(len(sp)); w=sp[8:40]; kk=k[8:40]; km=(w*kk).sum()/w.sum(); pk=8+np.argmax(w); half=w>0.5*w.max(); N=float(z['N'][0]); tot=w.sum(); o+='| %s power-weighted |k| (8-39) %.1f = %.1f cells, peak %d = %.1f cells, l %.0f, half-max bins %s'%(lab,km,N/km,pk,N/pk,18*pk,(int(kk[half][0]),int(kk[half][-1])))
        print(o)
