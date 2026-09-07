import sys,glob,numpy as np; sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/vis/python'); import bin_convert
N='/viper/u2/jinma/ATHENAK/bench/cs_mhd_prod_nan/bin/'
files=sorted(glob.glob(N+'*.bin')); prev=None
for f in files:
    r=bin_convert.read_binary(f); mb=r['mb_data']; g=np.asarray(r['mb_geometry'])
    rho=np.asarray(mb['dens']); e=np.asarray(mb['eint']); v=np.stack([np.asarray(mb[k]) for k in ('velx','vely','velz')]); b=np.stack([np.asarray(mb[k]) for k in ('bcc1','bcc2','bcc3')])
    bad=~(np.isfinite(rho)&np.isfinite(e)&np.isfinite(v).all(0)&np.isfinite(b).all(0))
    sp=np.sqrt((v**2).sum(0)); bb=np.sqrt((b**2).sum(0))
    line='%s rot %.3f: nan %d'%(f[-9:-4],r['time']/3.05e5,bad.sum())
    if not bad.any():
        m,k,j,i=np.unravel_index(np.nanargmax(sp),sp.shape); m2,k2,j2,i2=np.unravel_index(np.nanargmax(bb),bb.shape)
        line+=' | vmax %.1f km/s (blk %d i %d) Bmax %.0f G (blk %d i %d) rho min %.1e eint min %.2e'%(sp.max()/1e5,m,i,bb.max(),m2,i2,rho.min(),e.min())
        prev=(rho,e,sp,bb)
    else:
        m,k,j,i=np.where(bad); print(line)
        for mm in np.unique(m):
            sel=m==mm; print('   block %d (x2 %.2f..%.2f x3 %.2f..%.2f): %d nan; i %d..%d, j %d..%d, k %d..%d'%(mm,g[mm,2],g[mm,3],g[mm,4],g[mm,5],sel.sum(),i[sel].min(),i[sel].max(),j[sel].min(),j[sel].max(),k[sel].min(),k[sel].max()))
        if bad.sum()<2000:
            # print the nan cells and the previous dump's values there
            for idx in list(zip(m,k,j,i))[:40]:
                pr=prev
                print('     cell %s prev: rho %.2e eint %.2e |v| %.1f km/s |B| %.0f G'%(idx,pr[0][idx],pr[1][idx],pr[2][idx]/1e5,pr[3][idx]) if pr is not None else '     cell %s'%(idx,))
        break
    print(line)
