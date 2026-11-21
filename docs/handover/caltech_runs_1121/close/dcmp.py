import sys, numpy as np
# dcmp.py A.bin B.bin [icut]: compare VGD_DUMP_I files (5 int64 extents m,d,k,j,i, LayoutRight doubles)
def rd(f):
    sh=np.fromfile(f,dtype=np.int64,count=5); a=np.memmap(f,dtype=np.float64,mode='r',offset=40,shape=tuple(sh)); return a
a=rd(sys.argv[1]); b=rd(sys.argv[2]); ic=int(sys.argv[3]) if len(sys.argv)>3 else 0
print('shape',a.shape)
ni=a.shape[4]
for i in range(0,ni,max(1,ni//16)):
    x=np.asarray(a[...,i]); y=np.asarray(b[...,i])
    m=np.abs(x).max(); print(i,'maxA %.3e maxB %.3e nanA %d nanB %d maxrel %.2e'%(m,np.abs(y).max(),np.isnan(x).sum(),np.isnan(y).sum(),np.nanmax(np.abs(x-y))/max(m,1e-300)))
