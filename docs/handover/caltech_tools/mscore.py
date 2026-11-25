# mscore.py log.0.gz -> per-cycle SDOMOM ledger (poloidal momentum the floors removed) and the comparators
import gzip,re,sys
t=gzip.open(sys.argv[1],'rt',errors='replace').read()
tot={int(m[1]):[float(x) for x in m.groups()[1:]] for m in re.finditer(r'SDOMOM total cycle (-?\d+) dt (\S+) p_r (\S+) p_th (\S+) \|p_pol\| (\S+) KE_pol (\S+) M (\S+)',t)}
W={'update1':.5,'which2':.5,'update2':1.,'post2':1.,'flux3':.5,'other':1.,'flux3r':0.}   # end-of-step weight (flux3r: not persistent)
S={}
for m in re.finditer(r'SDOMOM sum cycle (-?\d+) (\S+) (\S+) cells (\d+) dp_r (\S+) \|dp_r\| (\S+) dp_th (\S+) \|dp_th\| (\S+) dKE (\S+) deint (\S+) dLz (\S+) \|p_pol\|pre (\S+) dE (\S+)',t):
    c,tag,cause=int(m[1]),m[2],m[3]; f=[float(x) for x in m.groups()[4:]]
    S.setdefault(cause,{}).setdefault(tag,[]).append((c,int(m[4]),f))
cyc=sorted(c for c in tot if c>=0); print('cycles',cyc[0],cyc[-1],len(cyc))
for cause in S:
    print('==',cause,'tag: cells/cycle | sum|dp_r| sum|dp_th| sum dp_r | dKE | max|dLz| max|dE| (raw; mean per cycle)')
    for tag,r in S[cause].items():
        n=len(cyc); a=lambda i:sum(x[2][i] for x in r)/n
        print('  %-8s %10.1f | %.3e %.3e %+.3e | %.3e | %.1e %.1e  %s'%(tag,sum(x[1] for x in r)/n,a(1),a(3),a(0),a(4),max(abs(x[2][6]) for x in r),max(abs(x[2][8]) for x in r),'NOT PERSISTENT' if tag=='flux3r' else 'x%.1f'%W[tag]))
L={}
for cause in S:
    L[cause]={c:sum(W[g]*(f[1]+f[3]) for g,r in S[cause].items() for cc,_,f in r if cc==c) for c in cyc}
    print('weighted |dp_r|+|dp_th| per cycle',cause,'mean %.3e max %.3e at cycle %d; share of |p_pol| mean %.2e max %.2e'%(sum(L[cause].values())/len(cyc),max(L[cause].values()),max(L[cause],key=L[cause].get),sum(L[cause][c]/tot[c][3] for c in cyc)/len(cyc),max(L[cause][c]/tot[c][3] for c in cyc)))
fc=[(int(c),float(p)) for c,p in re.findall(r'floorcap cycle (\d+) .*?\|dp\| (\S+)',t)]
if len(fc)>1: print('floorcap |dp| printed at',fc[0],fc[-1],'= %.3e per cycle'%((fc[-1][1]-fc[0][1])/max(fc[-1][0]-fc[0][0],1)))
m=re.search(r'dfloor_flux_limit: .*removed=(\S+?),',t)
if m: print('dfloor_flux_limit x1 momentum removed, this job: %s (/%d cycles = %.3e per cycle)'%(m[1],len(cyc),abs(float(m[1]))/len(cyc)))
else: print('dfloor_flux_limit line NOT FOUND')
p=[float(x) for x in re.findall(r'CLEAN FALLBACK cycle \d+:.*floor_pay (\S+)',t)]; print('floor_pay max %.3e'%max(p) if p else 'floor_pay NOT FOUND','CW',t.count('CLEAN WARNING'))
w=re.findall(r'SDOMOM worst cycle (\d+) tag (\d+) rank (\d+) .*? r (\S+) .*? dpV (\S+)',t); w=sorted(w,key=lambda x:-float(x[4])); print('largest single cells',w[:3])
