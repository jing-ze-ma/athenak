import sys, glob, numpy as np
sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/run'); import bin_convert
for d in sys.argv[1:]:
    fs=sorted(glob.glob(d+'/bin/*.bin'))
    print('==',d)
    for f in fs:
        try:
            raw=bin_convert.read_binary(f)
        except Exception as e:
            print(f,'ERR',e); continue
        mb=raw['mb_data']; t=raw['time']/ (2*np.pi/ (2*np.pi)) 
        vr=np.asarray(mb['velx'])[:,:,:,0]; rho=np.asarray(mb['dens'])[:,:,:,0]
        b=np.sqrt(sum(np.asarray(mb[k])[:,:,:,0]**2 for k in ['bcc1','bcc2','bcc3']))
        fin=np.isfinite(vr).all()
        n1=(vr<-1e5).sum(); n01=(vr<-1e4).sum()
        print(f'{f.split("/")[-1][-9:-4]} t={raw["time"]:.4e} finite={fin} n(vr<-1km/s)={n1} n(vr<-0.1)={n01} min vr={np.nanmin(vr)/1e5:.3f} km/s  Bmax(i=0)={np.nanmax(b)/np.sqrt(4*np.pi) if False else np.nanmax(b):.3e}  rho_min(i=0)={np.nanmin(rho):.3e}')
