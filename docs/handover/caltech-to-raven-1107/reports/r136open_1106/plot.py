import json,re,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
E='/resnick/groups/carnegie_poc/jingze/r136open_1106';P='/resnick/groups/carnegie_poc/jingze/r136a1_1104/runB/prod'
T0=432000.;D=86400.
def log(p):
    C=[]
    for l in open(p,errors='ignore'):
        m=re.match(r'elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)',l)
        if m: C.append([float(x) for x in m.groups()])
    return np.array(C)
col={'B':'#d62728','M':'#1f77b4','prod':'#555555','BL':'#2ca02c','ML':'#9467bd','C':'#8c564b'}
import os
ARMS=[a for a in ('B','M','BL','ML','C') if os.path.exists(E+'/an_%s.json'%a)]
fig,ax=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
a=ax[0,0]
for n,p in [('prod',P+'/run.log')]+[(a,E+'/runs/%s/run.log'%a) for a in ARMS]:
    c=log(p);c=c[(c[:,2]>=T0)&(c[:,2]<=T0+17280)]
    a.plot(c[:,2]/D,c[:,3],'-' if n!='prod' else 'o',ms=3,lw=.8,color=col[n],label={'prod':'control = production (same rst, 2 GPUs, every 50th cycle)','B':'B: no sponge, outflow','M':'M: sponge mode all (damps v_t)','BL':'BL: no sponge + linear source (diag. binary)','ML':'ML: mode all + linear source','C':'C: control noise member (E x (1+1e-10 u))'}[n])
a.set(xlabel='t [d]',ylabel='dt [s]',title='Time step per cycle',ylim=(0,7));a.legend(frameon=False,fontsize=9);a.grid(alpha=.3)
a=ax[0,1]
for n in ARMS:
    d=json.load(open(E+'/an_%s.json'%n));ks=sorted(d,key=int);t=[d[k]['t_d'] for k in ks]
    a.plot(t,[d[k]['50_90']['vr_rms'] for k in ks],'-o',color=col[n],label=n+': rms v_r 50-90 R$_\\odot$')
    a.plot(t,[-d[k]['top']['vr_med'] for k in ks],'--s',color=col[n],label=n+': -median v_r, top shell')
for n,f,k in (('b1 noise run (5.024 d)','b1','22'),('production (5.25 d)','prod','21')):
    q=json.load(open(E+'/an_%s.json'%f))[k];a.plot(q['t_d'],q['50_90']['vr_rms'],'k^' if f=='b1' else 'kv',label=n+': rms v_r 50-90')
a.set(xlabel='t [d]',ylabel='km/s',title='Atmosphere: infall without the sponge');a.legend(frameon=False,fontsize=8);a.grid(alpha=.3)
a=ax[1,0]
for n in ARMS:
    d=json.load(open(E+'/an_%s.json'%n));ks=sorted(d,key=int);t=[d[k]['t_d'] for k in ks]
    a.semilogy(t,[d[k]['35_50']['vtmax'] for k in ks],'-o',color=col[n],label=n+': max v_t, 35-50 R$_\\odot$')
    a.semilogy(t,[max(d[k]['ntop_T_27'],.5) for k in ks],':x',color=col[n],label=n+': cells at EOS table-top T (r>27)')
for f,k,m in (('b1','22','^'),('prod','21','v')):
    q=json.load(open(E+'/an_%s.json'%f))[k];a.semilogy(q['t_d'],q['35_50']['vtmax'],'k'+m,label={'b1':'b1 noise run','prod':'production'}[f]+': max v_t')
a.axhline(1e4,color='k',lw=.6,ls='--');a.text(5.002,1.1e4,'vceil 10$^4$ km/s',fontsize=8)
a.set(xlabel='t [d]',title='Hot fast cells (v_t is not damped by the radial sponge)');a.legend(frameon=False,fontsize=8);a.grid(alpha=.3)
a=ax[1,1];T=json.load(open(E+'/tcmp.json'));r=np.array(T['r'])
lab={'T0':'T0 sponge off (= B)','T1':'T1 + WB reconstruction everywhere','T2b':'T2b + isothermal walk','T3b':'T3b + no radiative kick','T5':'T5 + linear gravity source above 27.5 R$_\\odot$ (new key)','T6':'T6 new binary, no key'}
sty={'T0':('k','-'),'T1':('#2ca02c','--'),'T2b':('#9467bd',':'),'T3b':('#ff7f0e','-.'),'T5':('#1f77b4','-'),'T6':('#7f7f7f',':')}
for n in T['acc']:
    if n not in lab: continue
    v=np.array(T['acc'][n]['1080']);s=(r>36)
    a.plot(r[s],v[s],color=sty[n][0],ls=sty[n][1],lw=1.5,label=lab[n])
ph=np.array([(40,-1342),(45,-1079),(50,-871),(55,-717),(60,-623),(65,-547),(70,-489),(75,-427),(80,-385),(85,-349),(90,-307),(95,-283),(98,-268),(99.2,-260)])
a.plot(ph[:,0],ph[:,1],'r-',lw=2.5,alpha=.5,label='expected: -g(1-$\\Gamma$) (measured $\\Gamma$, $\\nabla P_{gas}$ ~ 0)')
a2=a.twinx();a2.semilogy(r[r>36],np.array(T['dr'])[r>36],color='#bbbbbb',lw=1);a2.set_ylabel('dr [R$_\\odot$] (grey)',color='#888888')
a.set(xlabel='r [R$_\\odot$]',ylabel='mean v_r / t after 1080 s [cm s$^{-2}$]',title='Force tests from rst 5.0 d, sponge off',xlim=(36,100),ylim=(-2500,0));a.legend(frameon=False,fontsize=8,loc='lower left');a.grid(alpha=.3)
fig.suptitle('R136a1 deep B, 5.00-5.20 d: sponge variants vs production control',fontsize=13)
fig.savefig(E+'/r136open.png',dpi=120);print('saved')
