import numpy as np, glob, json, sys
from PIL import Image
from scipy.ndimage import uniform_filter, binary_dilation
from multiprocessing import Pool
ST={'m35mms':140}
def ld(f):
    a=np.asarray(Image.open(f).convert('RGB').resize((960,540),Image.BILINEAR)).max(-1)>14; return a
def one(s):
    F=sorted(glob.glob('frames/%s_*.png'%s)); U=np.zeros((540,960),bool)
    from multiprocessing.pool import ThreadPool
    for m in ThreadPool(24).imap_unordered(ld,F): U|=m
    U[:60]=False; U[470:]=False; U[:,650:]=False; U[(540-115):,:240]=False           # title, v_r bar, T_eff bar, old scale bar (wedge is right of x = 240 px half-res there)
    U[(540-92):(540-70),340:500]=False           # the old bar of this layout (x 690-970, y 140-185)
    ys,xx=np.where(U); cx,cy=xx.mean(),ys.mean()
    bw=max(ST[s],150)//2+2; bh=30; g=13                                        # half-res box: bar or label width, label + bar height, gap 13 px = 26 px full = 2.4 % of the height
    occ=binary_dilation(U,structure=np.ones((2*g+1,2*g+1))); I=np.cumsum(np.cumsum(np.pad(occ,((1,0),(1,0))).astype(int),0),1); best=None
    for y in range(120,540-65-bh):                                               # y = top row (half-res) of the box
        for x in range(20,480):
            if x+bw/2>=cx or y+bh/2<=cy: continue
            if x<150 and y>480: continue
            if I[y+bh,x+bw]-I[y,x+bw]-I[y+bh,x]+I[y,x]>0: continue
            d=(x+bw/2-cx)**2+(y+bh/2-cy)**2
            if best is None or d<best[0]: best=(d,x,y)
    d,x,y=best; mid=2*(x+bw/2); yb=1080-2*(y+bh)+4                                # full-res bar: centre x, bar y from the bottom
    return s,(float(mid),float(yb))
if __name__=='__main__':
    with Pool(1) as p: R=dict(p.map(one,list(ST)))
    B=json.load(open('barpos.json')); B.update(R); json.dump(B,open('barpos.json','w')); print(R)
