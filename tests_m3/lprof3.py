import re, sys, math
fn=sys.argv[1]; which=sys.argv[2] if len(sys.argv)>2 else 'last'
blocks=[]; cur={}
for line in open(fn):
    m=re.match(r"rt_apply i=(\d+) .* Fb=([-\dEe.+]+) Ft=([-\dEe.+]+) .* Ab=([-\dEe.+]+) At=([-\dEe.+]+) ",line)
    if m:
        i=int(m.group(1))
        if i in cur: blocks.append(cur); cur={}
        cur[i]=tuple(float(m.group(g)) for g in (2,3,4,5))
if cur: blocks.append(cur)
F=blocks[0] if which=='first' else blocks[-1]
ks=sorted(F); rows=[]
for i in ks:
    fb,ft,ab,at=F[i]; rows.append((i,ab,fb,96.0*ab*fb))
i=ks[-1]; fb,ft,ab,at=F[i]; rows.append((i+1,at,ft,96.0*at*ft))
Ls=[x[3] for x in rows]
print("nblocks=%d  %s block: L(bot)=%.6e L(top)=%.6e  L max/min=%.6f  spread=%.4f %%"%(
  len(blocks),which,Ls[0],Ls[-1],max(Ls)/min(Ls),100.0*(max(Ls)-min(Ls))/(sum(Ls)/len(Ls))))
