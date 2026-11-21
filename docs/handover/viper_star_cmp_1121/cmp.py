import numpy as np,sys
runs={'m35mms (viper mesaR, general EOS)':('/viper/ptmp2/jinma/m35mms_1115/prod_viper_ours_mesaR/m35mms.user.hst',86400,[0,.05,.1,.2,.33,.5,.66,.8,1,1.5,2,3]),
 'm35mms (viper ours TOPS-R, general EOS)':('/viper/ptmp2/jinma/m35mms_1115/prod_viper_ours/m35mms.user.hst',86400,[0,.1,.2,.33,.5,.66,1]),
 't9l62 our-wall':('/viper/ptmp2/jinma/t9l62_1113/prod_viper/t9l62.user.hst',86400,[0,.25,.5,1,1.83,3,5,7.3,8.7]),
 't9l62 paper-BC':('/viper/ptmp2/jinma/t9l62_1113/prod_viper_paperbc/t9l62.user.hst',86400,[0,.5,1,1.83,3,5,7.9]),
 't19l64':('/viper/ptmp2/jinma/t19l64_1119/prod_viper/t19l64.user.hst',3600,[0,1,2,4,6,8,12,20,30,38,44])}
for k,(f,u,T) in runs.items():
    a=np.loadtxt(f); t,i=np.unique(a[:,0],return_index=True); a=a[i]; tt=a[:,0]/u
    print(f'== {k}: t_end {tt[-1]:.2f} {"d" if u==86400 else "h"}; col7 range {a[:,6].min():.4e}..{a[:,6].max():.4e}')
    Lin=np.median(a[:,6])
    for x in T:
        j=np.argmin(abs(tt-x))
        if abs(tt[j]-x)>0.2*max(x,0.05): continue
        print(f'  t {tt[j]:6.2f}  dt {a[j,1]:7.2f}  KE {a[j,10]:9.2e}  Ltop/Lin {a[j,5]/Lin:6.3f}  Lbot/Lin {a[j,2]/Lin:6.3f}')
