import sys,numpy as np
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/doface_1105/page/pipe')
import bin_convert_repro as bc, bsg_figs as B
D='/resnick/groups/carnegie_poc/jingze/doface_1105/prod_vsnmlat_ns/'
d=bc.read_binary(D+'bin/bsg3d.hydro_w.00034.bin'); a=B.assemble(d)
print(d['var_names'])
g=np.loadtxt(D+'bsg3d.x1grid.txt',comments='#',max_rows=262); x1v=g[3:259,2]; dx1=g[3:259,3]
nx3,nx2,nx1=a['dens'].shape
th=np.linspace(np.pi/3,2*np.pi/3,nx2+1); thv=0.5*(th[1:]+th[:-1]); dth=th[1]-th[0]; dph=(np.pi/3)/nx3
rho=a['dens'].astype(np.float64)
P=a['pres'] if 'pres' in a else (2/3)*a['eint']
cs=np.sqrt(5/3*P.astype(np.float64)/rho)
v1,v2,v3=[np.abs(a[k].astype(np.float64)) for k in ('velx','vely','velz')]
R=x1v[None,None,:]; S=np.sin(thv)[None,:,None]
rsun=6.957e10; merge=np.where(x1v<2.43495e12,2.0,1.0)[None,None,:]
for lab,mf in (('unmerged',1.0),('merged x2 below 35',merge)):
    t1=dx1[None,None,:]/(v1+cs); t2=mf*R*dth/(v2+cs); t3=mf*R*S*dph/(v3+cs)
    for nm,t in (('x1',t1),('x2',t2),('x3',t3)):
        i=np.unravel_index(np.argmin(t),t.shape)
        print(f"{lab:20s} {nm}: cfl0.3*min dt = {0.3*t[i]:7.1f} s at r={x1v[i[2]]/rsun:5.2f} Rsun  cs={cs[i]/1e5:6.1f} km/s |v|={[v1,v2,v3][int(nm[1])-1][i]/1e5:6.1f}")
    tm=np.minimum(np.minimum(t1,t2),t3)*0.3
    prof=tm.min(axis=(0,1)); 
    for r0 in (20.5,25,30,34.5,35.5,40,45,50,55,60,70,79):
        j=np.argmin(abs(x1v/rsun-r0)); print(f"   r={x1v[j]/rsun:5.1f}: min dt {prof[j]:7.1f} s")
