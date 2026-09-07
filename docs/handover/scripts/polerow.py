import sys, numpy as np
sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/run'); import bin_convert
for f in sys.argv[1:]:
    raw=bin_convert.read_binary(f); mb=raw['mb_data']; t=raw['time']; g=np.asarray(raw['mb_geometry'])
    B=[np.asarray(mb[k]) for k in ['bcc1','bcc2','bcc3']]; rho=np.asarray(mb['dens']); vr=np.asarray(mb['velx'])
    north=[m for m in range(g.shape[0]) if g[m,2]<1e-6]; south=[m for m in range(g.shape[0]) if abs(g[m,3]-np.pi)<1e-6]
    nj=B[0].shape[2]
    print(f'{f}: rot={t/3.05e5:.2f}  north MBs {north} south {south}')
    for name,mbs,rows in [('N',north,[0,1,2,3,5,8]),('S',south,[nj-1,nj-2,nj-3,nj-4,nj-6,nj-9])]:
        for i in [0,1,4,16,48]:
            s=[]
            for j in rows:
                br=np.concatenate([B[0][m,:,j,i] for m in mbs]); bt=np.concatenate([B[1][m,:,j,i] for m in mbs]); bp=np.concatenate([B[2][m,:,j,i] for m in mbs])
                s.append(f'j{j}:{np.abs(br).max():.0f}/{np.abs(bt).max():.0f}/{np.abs(bp).max():.0f}')
            print(f'  {name} i={i:2d} max|Br|/|Bt|/|Bp|  '+'  '.join(s))
        # spectrum in phi of Br at the polar row, i=0 and i=4
        for i in [0,4]:
            j=rows[0]; br=np.concatenate([B[0][m,:,j,i] for m in mbs]); n=len(br); F=np.abs(np.fft.rfft(br-br.mean()))**2; tot=F.sum()+1e-300
            top=np.argsort(F)[::-1][:3]
            print(f'  {name} i={i} j={j} Br pole row: mean {br.mean():.1f} rms {br.std():.1f}; power frac m=1 {F[1]/tot:.2f} m=2 {F[2]/tot:.2f} Nyq(m={n//2}) {F[n//2]/tot:.2f}; top m {top.tolist()} rho min {np.concatenate([rho[m,:,j,i] for m in mbs]).min():.2e} vr min/max {np.concatenate([vr[m,:,j,i] for m in mbs]).min()/1e5:.2f}/{np.concatenate([vr[m,:,j,i] for m in mbs]).max()/1e5:.2f} km/s')
