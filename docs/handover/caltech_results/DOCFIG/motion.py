# frame-to-frame change of the mapped fields: correlation time of v_r and T_eff in the maps, per dump pair (1,2,4 dumps), and peak shift [deg]
import sys; sys.argv=['x','bsg','none']
import numpy as np, docfig as d
for s in ('m35mms','m35zams'):
    F=d.files(s); out=[]
    for a in (150,250,350,450):
        for k in (1,2,4):
            A=d.load(s,F[a]); B=d.load(s,F[a+k]); r=[]
            for q in ('vr','te','sl'):
                x=np.ma.masked_invalid(A[q]).filled(0).astype(float); y=np.ma.masked_invalid(B[q]).filled(0).astype(float)
                x=x-x.mean(); y=y-y.mean(); r.append(float((x*y).sum()/np.sqrt((x*x).sum()*(y*y).sum())))
            out.append((a,k,[round(v,3) for v in r]))
    print(s,'dt per dump [min] %.1f'%((d.load(s,F[451])['t']-d.load(s,F[450])['t'])*1440)); [print('  dump',o[0],'+%d'%o[1],'corr v_r, T_eff, drho slice',o[2]) for o in out]
