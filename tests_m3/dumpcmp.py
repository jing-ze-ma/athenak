import re,sys
fn=sys.argv[1]
rows=[]
for l in open(fn):
    m=re.search(r"rt_col3_flux i=(\d+) dtau=([-\dEe.+]+) w=[-\dEe.+]+ srcdx=([-\dEe.+]+) divF=([-\dEe.+]+) dif=([-\dEe.+]+) sw_srcdx=([-\dEe.+]+) F3lo/A=([-\dEe.+]+) Fblo=([-\dEe.+]+)",l)
    if m: rows.append([float(x) for x in m.groups()])
if not rows: print("no dump"); sys.exit()
def rel(a,b): 
    s=max(abs(a),abs(b));  return abs(a-b)/s if s>0 else 0.0
msrc=max(rel(r[2],r[5]) for r in rows)
mtel=max(rel(r[2],r[3]) for r in rows)
mflx=max(rel(r[6],r[7]) for r in rows)
print("cells=%d  tau/cell %.3g..%.3g" % (len(rows), rows[0][1], rows[-1][1]))
print("  max rel |srcdx - sw_srcdx|  (mode-3 source vs the sweep's) = %.3e" % msrc)
print("  max rel |srcdx - divF|      (flux telescoping)             = %.3e" % mtel)
print("  max rel |F3lo/A - Fblo|     (mode-3 face flux vs sweep's)  = %.3e" % mflx)
