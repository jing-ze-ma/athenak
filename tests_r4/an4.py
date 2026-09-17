import sys, numpy as np
R = 1.18585e11/0.50           # 2.3717e11 cm, the photosphere; x1min = 0.50 R
print("%-10s %6s %8s %8s %8s %8s"%("dir","r/R","F2s/Freq","Frd/Freq","ratio","w"))
for p in sys.argv[1:]:
    d=np.loadtxt(p+"/mltfaces.txt")
    r=d[:,1]/R; freq,frd,f2s,w=d[:,10],d[:,11],d[:,17],d[:,12]
    sel=(r>=0.55)&(r<=0.97)&(freq!=0)
    a=f2s[sel]/freq[sel]; b=frd[sel]/freq[sel]
    rat=a/b
    print("### %s  N=%d"%(p,sel.sum()))
    print("   F2s/Freq  min %.4f max %.4f median %.4f"%(a.min(),a.max(),np.median(a)))
    print("   Frd/Freq  min %.4f max %.4f median %.4f"%(b.min(),b.max(),np.median(b)))
    print("   F2s/Frd   min %.4f max %.4f median %.4f  max|1-x| %.4f  wmax %.2e"
          %(rat.min(),rat.max(),np.median(rat),np.abs(rat-1).max(),w.max()))
    rr=r[sel]
    for k in range(0,sel.sum(),4):
        print("   %6.3f %8.4f %8.4f %8.4f"%(rr[k],a[k],b[k],rat[k]))
