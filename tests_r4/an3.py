import sys, numpy as np
def rin_of(tag):
    return 1.0e12
rows=[]
for tag in ["n0_t100_r2.0","n1_t100_r2.0","n3_t100_r2.0","n7_t100_r2.0",
            "n3_t30_r2.0","n3_t300_r2.0","n3_t100_r1.1","n3_t100_r1.5"]:
    d=np.loadtxt(sys.argv[1]+"/"+tag+"/mltfaces.txt")
    r,T=d[:,1],d[:,2]; freq,fdiff,f2s,w=d[:,10],d[:,11],d[:,17],d[:,12]
    # cell optical depth from the local kappa*rho*dr is not dumped; use dtau proxy via
    # the analytic requirement: restrict to the interior faces 10..118 as the README did
    sel=np.arange(6,len(r)-8)
    ratio=f2s[sel]/freq[sel]; rd=fdiff[sel]/freq[sel]
    # keep only faces where the Rosseland control itself is within 1% (dtau_cell > 0.3)
    ok=np.abs(rd-1.0)<0.01
    rr=ratio[ok]
    rows.append((tag,len(rr),rr.min(),rr.max(),np.median(rr),
                 np.abs(rr-1).max(), w.max()))
print("%-16s %5s %8s %8s %8s %10s %8s"%("case","N","min","max","median","max|1-x|","wmax"))
for t in rows: print("%-16s %5d %8.4f %8.4f %8.4f %10.4f %8.1e"%t)
