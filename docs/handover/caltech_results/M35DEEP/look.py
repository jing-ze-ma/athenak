import glob, numpy as np, sys
for run in sys.argv[1:]:
    S=[dict(np.load(f)) for f in sorted(glob.glob('stars/%s/s_*.npz'%run))]; S.sort(key=lambda d: float(d['t'])); z=S[-1]; r=z['r']; N3,N2=z['N']; ex=z['ext']
    print('=====',run,'dumps',len(S),'t %.3e-%.3e'%(S[0]['t'],z['t']),'N3 x N2',N3,N2,'widths deg',np.degrees(ex[1]-ex[0]),np.degrees(ex[3]-ex[2]),'dtype',z['dtype'])
    P=np.load('stars/prof_%s.npz'%run) if run in ('m35zams','bsg','t9l62','t9mlt','m35mms') else None
    ni=len(r)
    for i in range(1,ni,max(1,ni//16)):
        sp=z['sp'][:,i]; sa=z['sa'][:,i]; tot=sp[1:].sum(); bands=[(1,3),(4,8),(9,15),(16,19),(20,28),(29,40),(41,64),(65,128)]
        k=np.arange(len(sp)); pk=3+np.argmax(sp[3:100]); tp=z['top'][i]
        print(' r %.4e (%.3f) rms %.2e cm/s | band shares %%:'%(r[i],r[i]/r[0],z['rms'][i]),' '.join('%d-%d:%.1f'%(a,b,100*sp[a:b+1].sum()/tot) for a,b in bands),'| axis share of 16-40: %.2f (count %.2f)'%(sa[16:41].sum()/max(sp[16:41].sum(),1e-99), z['cax'][16:41].sum()/z['cnt'][16:41].sum()),'| peak bin (>=3) %d = %.1f cells | top modes (kx,ky,k256,share):'%(pk,256/pk),[(int(a),int(b),int(c),round(float(s),2)) for a,b,c,s in tp[:3]], ('| N %.2e'%np.sqrt(max(np.interp(r[i],P['r'],P['N2']),0)) if P is not None else ''), '| vt/vr %.2f'%(z['vt'][i]/max(z['rms'][i],1e-30)))
