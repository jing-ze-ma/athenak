import re,sys
rows=[]
for l in open(sys.argv[1]):
    m=re.search(r"rt_col3_flux i=(\d+) dtau=([-\dEe.+]+) w=[-\dEe.+]+ srcdx=([-\dEe.+]+) divF=([-\dEe.+]+) dif=([-\dEe.+]+) sw_srcdx=([-\dEe.+]+) F3lo/A=([-\dEe.+]+) Fblo=([-\dEe.+]+)",l)
    if m: rows.append([float(x) for x in m.groups()])
segs=[];cur=[];prev=None
for v in rows:
    if prev is not None and v[0]<=prev: segs.append(cur);cur=[]
    cur.append(v);prev=v[0]
if cur: segs.append(cur)
s=segs[0]; S=max(abs(r[2]) for r in s)
print("  %d cells, tau/cell %.3g (bottom) .. %.3g (top)"%(len(s),s[0][1],s[-1][1]))
print("  max |srcdx - sw_srcdx| / max|srcdx| = %.3e"%(max(abs(r[2]-r[5]) for r in s)/S))
print("  max |srcdx - divF|     / max|srcdx| = %.3e"%(max(abs(r[2]-r[3]) for r in s)/S))
