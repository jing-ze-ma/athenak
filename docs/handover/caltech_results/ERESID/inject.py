import sys, numpy as np; sys.argv=['x']
src=open('/resnick/groups/carnegie_poc/jingze/eres_1124/eres.py').read().split("P = []; pr =")[0]; exec(src)
rng=np.random.default_rng(1)
for star,d,base,sw_ in STARS:
    ch=rhst(J+d+'/'+base+'.user.hst'); n,a=ch[0]; g=lambda nm:a[:,n.index(nm)]
    t=g('time'); E=g('Etot'); M=g('M_tot'); rate=g('L_in')-g('L_top')+g('P_esrc')-g('P_sponge')-g('E_hflx'); dt=np.diff(t)
    R=(np.diff(E)-dt*rate[1:])/E[1:]; RM=(np.diff(M)+dt*g('M_hflx')[1:])/M[1:]; res=R-rmed(R); mad=np.median(abs(res)); rm=RM-rmed(RM); madm=np.median(abs(rm))
    print('%s: rows %d MAD E %.2e M %.2e'%(star,len(R),mad,madm), end=' | ')
    for J_ in (1e-9,1e-10,1e-12,1e-13,1e-14):
        k=rng.choice(len(R)-200,20,replace=False)+100; x=R.copy(); x[k]+=J_*rng.choice((-1,1),20); z=x-rmed(x); m=np.median(abs(z)); hit=(abs(z[k])>10*m).sum()
        xm=RM.copy(); xm[k]+=J_*rng.choice((-1,1),20); zm=xm-rmed(xm); mm=np.median(abs(zm)); hm=(abs(zm[k])>10*mm).sum()
        print('jump %.0e: E flagged %d/20, M %d/20;'%(J_,hit,hm), end=' ')
    print()
