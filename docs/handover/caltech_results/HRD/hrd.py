"""hrd.py: HR-diagram table. L_top and L_in of the user hst (wedge) x 4 pi/Omega; T_eff = median over columns of the comoving-flux T_eff at tau_R = 2/3 of each column
(docfig_1125/p23.py, the page definition); R_ph = median over columns of the tau_R = 2/3 radius; T_eff,LR = (L_top_sphere/(4 pi R_ph^2 sigma))^1/4. Read only."""
import numpy as np, glob, csv
G='/resnick/groups/carnegie_poc/jingze/'; LS=3.828e33; RS=6.957e10; SIG=5.670374419e-5
S=dict(bsg=('prod_c1_1124/bsg/bsg3d',30,57.01,np.pi/3,np.pi/3,2*np.pi/3,5.161),m35mms=('prod_c1_1124/m35mms/m35mms',5,10.01,0.34906585039886590,1.39626340159546358,1.74532925199432953,5.444),m35zams=('prod_c1_1124/m35zams/m35zams',5,10.01,0.064,1.5387520817282807,1.6028405718615124,5.238),t9mlt=('prod_mlt_1125/t9mlt/t9l62',None,None,np.pi,0.78539816339744831,2.3561944901923448,5.811))
def st(x): return dict(mean=float(np.mean(x)),sd=float(np.std(x)),lo=float(np.min(x)),hi=float(np.max(x)))
rows=[]
for s,(pre,t0,t1,dp,a,b,dl) in S.items():
    om=dp*(np.cos(a)-np.cos(b)); ws=4*np.pi/om
    hdr=[l for l in open(G+pre+'.user.hst') if '[1]=time' in l][0]; cols=[c.split('=')[1] for c in hdr[1:].split()]; H=np.array([list(map(float,l.split())) for l in open(G+pre+'.user.hst') if not l.startswith('#') and len(l.split())==len(cols)]); t=H[:,0]/86400; Lt=H[:,cols.index('L_top')]*ws/LS; Li=H[:,cols.index('L_in')]*ws/LS
    F=sorted(glob.glob(G+'docfig_1125/p23/%s_0*.npz'%s)); D=[]
    for f in F:
        z=np.load(f); D.append((float(z['t']),float(np.nanmedian(z['te'])),float(np.nanmedian(z['r']))/RS,int(z['nnone']),float(np.mean(~np.isfinite(z['te'])))))
    D=np.array(D); tl=D[-1,0]; t0=tl-5 if t0 is None else t0; t1=tl+.001 if t1 is None else t1; m=(D[:,0]>=t0-1e-6)&(D[:,0]<=t1); d=D[m]; h=(t>=t0)&(t<=t1)
    Tlr=(np.interp(d[:,0],t,Lt)*LS/(4*np.pi*(d[:,2]*RS)**2*SIG))**.25
    r=dict(star=s,t0=t0,t1=float(d[-1,0]),ndumps=int(m.sum()),last_hst=float(t[-1]),Omega=om,fourpi_over_Omega=ws,logL_design=dl,logL_imposed=float(np.log10(np.mean(Li[h]))),**{'L_'+k:v for k,v in st(Lt[h]).items()},**{'T_'+k:v for k,v in st(d[:,1]).items()},**{'TLR_'+k:v for k,v in st(Tlr).items()},**{'R_'+k:v for k,v in st(d[:,2]).items()},nan_frac_max=float(d[:,4].max()),ratio_T=float(np.mean(d[:,1]/Tlr)))
    r['logL_mean']=float(np.log10(r['L_mean'])); r['logL_sd']=float(np.std(np.log10(Lt[h]))); r['logL_lo']=float(np.log10(Lt[h].min())); r['logL_hi']=float(np.log10(Lt[h].max())); r['logR_mean']=float(np.log10(r['R_mean'])); r['logT_mean']=float(np.log10(r['T_mean'])); n=len(d)//2; r['T_trend']=float(d[n:,1].mean()-d[:n,1].mean()); r['R_trend']=float(d[n:,2].mean()-d[:n,2].mean()); r['logL_hrow_first']=float(np.log10(np.interp(t0,t,Lt))); rows.append(r)
    print(s,{k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()})
K=list(rows[0]); w=csv.DictWriter(open('hrd_table.csv','w'),K); w.writeheader(); w.writerows(rows)
