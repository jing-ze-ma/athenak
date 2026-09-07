"""Cubed-sphere dhj analysis: EOS inversion (T, p), panel -> lat/lon, zonal wind, regridding."""
import sys, numpy as np
sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/run'); import bin_convert

# ---------------- EOS table ----------------
class EOS:
    def __init__(self, fn):
        with open(fn) as f:
            lines=[next(f) for _ in range(6)]
        nx,ny,xmin,dx,ymin,dy=[float(v) for v in lines[2][1:].split()]
        self.nx,self.ny=int(nx),int(ny); self.xmin,self.dx,self.ymin,self.dy=xmin,dx,ymin,dy
        d=np.loadtxt(fn)
        self.le=d[:,0].reshape(self.ny,self.nx); self.lp=d[:,1].reshape(self.ny,self.nx)
        self.lt=self.ymin+self.dy*np.arange(self.ny); self.ld=self.xmin+self.dx*np.arange(self.nx)
    def invert(self, rho, eint):
        """rho, eint (erg/cm3) arrays -> T [K], p [erg/cm3]. Bilinear in log rho, linear-interp in T."""
        rho=np.asarray(rho,float); eint=np.asarray(eint,float); sh=rho.shape; rho=rho.ravel(); eint=eint.ravel()
        x=(np.log10(rho)-self.xmin)/self.dx; ix=np.clip(np.floor(x).astype(int),0,self.nx-2); fx=np.clip(x-ix,0,1)
        # e/rho curves at the two bracketing densities, interpolated in x -> one curve per cell (ny)
        le=self.le[:,ix]*(1-fx)+self.le[:,ix+1]*fx      # (ny, N)
        lp=self.lp[:,ix]*(1-fx)+self.lp[:,ix+1]*fx
        target=np.log10(eint/rho)
        fin=np.isfinite(le)
        lef=np.where(fin, le, -np.inf)
        # first finite row per column
        first=np.argmax(fin,axis=0)
        cnt=(lef<target[None,:]).sum(axis=0)      # rows below target (NaN rows count as below via -inf) 
        # subtract the NaN rows that were counted
        nnan=(~fin).sum(axis=0)
        k=cnt  # index of first row >= target ; rows [0,first) are NaN and counted in cnt (since -inf<target)
        k=np.clip(k, first+1, self.ny-1)
        e0=le[k-1,np.arange(len(k))]; e1=le[k,np.arange(len(k))]
        f=np.clip((target-e0)/(e1-e0),0,1)
        lT=self.lt[k-1]+f*self.dy
        lP=lp[k-1,np.arange(len(k))]*(1-f)+lp[k,np.arange(len(k))]*f
        T=10**lT; p=10**lP*rho
        return T.reshape(sh), p.reshape(sh)

# ---------------- cubed sphere geometry ----------------
FRAMES={0:([0,1,0],[0,0,1],[1,0,0]),1:([-1,0,0],[0,0,1],[0,1,0]),2:([0,-1,0],[0,0,1],[-1,0,0]),
        3:([0,1,0],[-1,0,0],[0,0,1]),4:([1,0,0],[0,0,1],[0,-1,0]),5:([0,1,0],[1,0,0],[0,0,-1])}
def panel_map(p, xi, eta):
    a,b,n=[np.array(v,float) for v in FRAMES[p]]
    x=np.tan(xi); y=np.tan(eta); dl=np.sqrt(1+x*x+y*y)
    X=(a[None,None,:]*x[...,None]+b[None,None,:]*y[...,None]+n[None,None,:])/dl[...,None]
    return X
def cs_geometry(raw):
    """returns per-MB arrays (nmb, nk, nj): lat, lon_ss (from substellar, -x), unit vectors exi, eeta, and zonal/meridional projections"""
    g=np.asarray(raw['mb_geometry']); nmb=g.shape[0]; nj=raw['nx2_mb']; nk=raw['nx3_mb']
    npan=6; per=nmb//npan
    out={k:np.zeros((nmb,nk,nj)) for k in ['lat','lon','exi_z','exi_m','eeta_z','eeta_m','cos']}
    for m in range(nmb):
        p=m//per
        xi=np.pi/4*(g[m,2]+(g[m,3]-g[m,2])*(np.arange(nj)+0.5)/nj)
        eta=np.pi/4*(g[m,4]+(g[m,5]-g[m,4])*(np.arange(nk)+0.5)/nk)
        XI,ETA=np.meshgrid(xi,eta)   # (nk,nj)
        X=panel_map(p,XI,ETA); h=1e-6
        exi=panel_map(p,XI+h,ETA)-panel_map(p,XI-h,ETA); exi/=np.linalg.norm(exi,axis=-1)[...,None]
        eeta=panel_map(p,XI,ETA+h)-panel_map(p,XI,ETA-h); eeta/=np.linalg.norm(eeta,axis=-1)[...,None]
        lat=np.arcsin(X[...,2]); lon=np.arctan2(-X[...,1],-X[...,0])   # 0 at substellar (-x)
        zon=np.stack([-X[...,1],X[...,0],np.zeros_like(lat)],-1); zon/=np.maximum(np.linalg.norm(zon,axis=-1),1e-30)[...,None]
        mer=np.stack([-X[...,2]*np.cos(np.arctan2(X[...,1],X[...,0])),-X[...,2]*np.sin(np.arctan2(X[...,1],X[...,0])),np.cos(lat)],-1)  # northward
        out['lat'][m]=lat; out['lon'][m]=lon; out['cos'][m]=(exi*eeta).sum(-1)
        out['exi_z'][m]=(exi*zon).sum(-1); out['eeta_z'][m]=(eeta*zon).sum(-1)
        out['exi_m'][m]=(exi*mer).sum(-1); out['eeta_m'][m]=(eeta*mer).sum(-1)
    return out
def winds(raw, geo):
    mb=raw['mb_data']; v2=np.asarray(mb['vely']); v3=np.asarray(mb['velz'])
    u=v2*geo['exi_z'][...,None]+v3*geo['eeta_z'][...,None]
    v=v2*geo['exi_m'][...,None]+v3*geo['eeta_m'][...,None]
    return u,v
def bfield(raw, geo):
    """bcc: (B_r, B.exi, s*B^eta) orthonormal (exi, eperp). eperp = (eeta - c exi)/s."""
    mb=raw['mb_data']; b1=np.asarray(mb['bcc1']); b2=np.asarray(mb['bcc2']); b3=np.asarray(mb['bcc3'])
    c=geo['cos'][...,None]; s=np.sqrt(1-c*c)
    # eperp . zon = (eeta.zon - c exi.zon)/s
    ez=(geo['eeta_z'][...,None]-c*geo['exi_z'][...,None])/s; em=(geo['eeta_m'][...,None]-c*geo['exi_m'][...,None])/s
    bz=b2*geo['exi_z'][...,None]+b3*ez; bm=b2*geo['exi_m'][...,None]+b3*em
    return b1,bz,bm,np.sqrt(b1*b1+b2*b2+b3*b3)

def level_interp(q, logp, levels):
    """q, logp: (nmb,nk,nj,ni) with p decreasing outward; levels: log10 p targets -> (nlev, nmb,nk,nj)"""
    out=np.full((len(levels),)+q.shape[:3], np.nan)
    lp=-logp  # increasing with i
    for l,L in enumerate(levels):
        t=-L
        idx=(lp<t).sum(axis=-1)  # number of cells with p > level
        ok=(idx>0)&(idx<q.shape[-1])
        i1=np.clip(idx,1,q.shape[-1]-1); i0=i1-1
        a=np.take_along_axis(lp,i0[...,None],-1)[...,0]; b=np.take_along_axis(lp,i1[...,None],-1)[...,0]
        qa=np.take_along_axis(q,i0[...,None],-1)[...,0]; qb=np.take_along_axis(q,i1[...,None],-1)[...,0]
        f=np.clip(np.where(b>a,(t-a)/np.where(b>a,b-a,1.0),0.5),0.0,1.0); val=qa+f*(qb-qa); out[l]=np.where(ok&(b>a),val,np.nan)
    return out
def latlon_bin(field, lat, lon, dlat=6, dlon=6):
    """field (..., nmb,nk,nj) -> mean on a lat/lon grid; returns grid centres and array (..., nlat, nlon)"""
    lats=np.arange(-90,90,dlat); lons=np.arange(-180,180,dlon)
    il=np.clip(((np.degrees(lat)+90)//dlat).astype(int),0,len(lats)-1); io=np.clip(((np.degrees(lon)+180)//dlon).astype(int),0,len(lons)-1)
    flat=field.reshape(field.shape[:-3]+(-1,)); il=il.ravel(); io=io.ravel()
    key=il*len(lons)+io; n=len(lats)*len(lons)
    cnt=np.bincount(key,minlength=n).astype(float)
    res=np.zeros(field.shape[:-3]+(n,))
    for idx in np.ndindex(field.shape[:-3]):
        f=flat[idx]; good=np.isfinite(f)
        res[idx]=np.bincount(key[good],weights=f[good],minlength=n)/np.maximum(np.bincount(key[good],minlength=n),1)
        empty=np.bincount(key[good],minlength=n)==0
        if empty.any() and (~empty).any():
            # nearest-neighbour fill of empty bins (polar rows, sp stretched rows) on the sphere
            LA,LO=np.meshgrid(np.radians(lats+dlat/2),np.radians(lons+dlon/2),indexing='ij'); LA=LA.ravel(); LO=LO.ravel()
            cx=np.cos(LA)*np.cos(LO); cy=np.cos(LA)*np.sin(LO); cz=np.sin(LA)
            gi=np.where(~empty)[0]; ei=np.where(empty)[0]
            dots=cx[ei,None]*cx[None,gi]+cy[ei,None]*cy[None,gi]+cz[ei,None]*cz[None,gi]
            res[idx][ei]=res[idx][gi[np.argmax(dots,axis=1)]]
    return lats+dlat/2, lons+dlon/2, res.reshape(field.shape[:-3]+(len(lats),len(lons)))
def zonal_mean(field, lat, dlat=6):
    lats=np.arange(-90,90,dlat); il=np.clip(((np.degrees(lat)+90)//dlat).astype(int),0,len(lats)-1).ravel()
    flat=field.reshape(field.shape[:-3]+(-1,)); res=np.zeros(field.shape[:-3]+(len(lats),))
    for idx in np.ndindex(field.shape[:-3]):
        f=flat[idx]; good=np.isfinite(f)
        c=np.bincount(il[good],minlength=len(lats)); res[idx]=np.bincount(il[good],weights=f[good],minlength=len(lats))/np.maximum(c,1); res[idx][c==0]=np.nan
    return lats+dlat/2, res

# ---------------- spherical polar geometry ----------------
def sp_geometry(raw, fstr=3.0):
    g=np.asarray(raw['mb_geometry']); nmb=g.shape[0]; nj=raw['nx2_mb']; nk=raw['nx3_mb']
    lat=np.zeros((nmb,nk,nj)); lon=np.zeros((nmb,nk,nj))
    for m in range(nmb):
        th=g[m,2]+(g[m,3]-g[m,2])*(np.arange(nj)+0.5)/nj
        xi=th/np.pi; th=np.pi/2*(1+np.sinh(fstr*(2*xi-1))/np.sinh(fstr))
        ph=g[m,4]+(g[m,5]-g[m,4])*(np.arange(nk)+0.5)/nk
        TH,PH=np.meshgrid(th,ph); lat[m]=np.pi/2-TH
        x=np.sin(TH)*np.cos(PH); y=np.sin(TH)*np.sin(PH); lon[m]=np.arctan2(-y,-x)
    return {'lat':lat,'lon':lon}
def sp_winds(raw):
    mb=raw['mb_data']; return np.asarray(mb['velz']), -np.asarray(mb['vely'])   # zonal, northward

def load(path, eos, kind):
    raw=bin_convert.read_binary(path); mb=raw['mb_data']
    T,p=eos.invert(np.asarray(mb['dens']),np.asarray(mb['eint']))
    if kind=='cs':
        geo=cs_geometry(raw); u,v=winds(raw,geo)
    else:
        geo=sp_geometry(raw); u,v=sp_winds(raw)
    d={'raw':raw,'rot':raw['time']/3.05e5,'T':T,'p':p,'lp':np.log10(p),'u':u,'v':v,'lat':geo['lat'],'lon':geo['lon'],'geo':geo,'rho':np.asarray(mb['dens'])}
    if 'bcc1' in mb:
        if kind=='cs': d['Br'],d['Bz'],d['Bm'],d['Babs']=bfield(raw,geo)
        else:
            d['Br']=np.asarray(mb['bcc1']); d['Bz']=np.asarray(mb['bcc3']); d['Bm']=-np.asarray(mb['bcc2']); d['Babs']=np.sqrt(d['Br']**2+d['Bz']**2+d['Bm']**2)
    return d
