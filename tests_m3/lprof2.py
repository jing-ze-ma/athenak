import re, sys, math
fn = sys.argv[1]
F = {}
for line in open(fn):
    m = re.match(r"rt_apply i=(\d+) .* Fb=([-\dEe.+]+) Ft=([-\dEe.+]+) .* Ab=([-\dEe.+]+) At=([-\dEe.+]+) ", line)
    if m:
        i = int(m.group(1))
        if i not in F: F[i] = tuple(float(m.group(g)) for g in (2,3,4,5))
ks = sorted(F)
rows = []
for i in ks:
    fb,ft,ab,at = F[i]
    rows.append((i, ab, fb, 96.0*ab*fb))
i=ks[-1]; fb,ft,ab,at = F[i]
rows.append((i+1, at, ft, 96.0*at*ft))
print("#    i      A[cm^2]       F[erg/cm^2/s]     L=A F * 96")
for k,(i,a,f,L) in enumerate(rows):
    if k % 8 == 0 or k == len(rows)-1:
        print("%6d  %.5e  %.8e  %.8e" % (i,a,f,L))
Ls=[x[3] for x in rows]
pos=[x for x in Ls if x>0]
print("# L max/min = %.5f ; spread/mean = %.4f %%" % (max(Ls)/min(Ls) if min(Ls)!=0 else float('nan'),
      100.0*(max(Ls)-min(Ls))/(sum(Ls)/len(Ls))))
