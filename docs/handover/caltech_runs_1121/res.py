import sys,glob,os,re
D,M=sys.argv[1:3]
L=open(D+'/a/run.log',errors='ignore').read().splitlines()
el=[(float(m.group(1)),int(m.group(2)),float(m.group(3))) for l in L for m in [re.search(r'elapsed=(\S+) cycle=(\d+) time=\S+ dt=(\S+)',l)] if m]
if len(el)>2:
    e=[x[0] for x in el]; d=[e[i+1]-e[i] for i in range(len(e)-1)]
    st=d[1:-1] if len(d)>2 else d
    print('cycles',el[0][1],'-',el[-1][1],'dt',sorted(set(x[2] for x in el)),'s/cycle(all)',['%.2f'%x for x in d],'steady(excl first+last) %.3f'%(sum(st)/len(st)),'first elapsed %.1f'%e[0],'last elapsed %.1f'%e[-1])
print([l for l in L if 'main loop time' in l][-1:]); print('CW a',sum('CLEAN WARNING' in l for l in L),'FATAL',sum('FATAL' in l for l in L), 'retries',[l for l in L if 'retries' in l][-1:] )
mem={}
for l in open(D+'/mem.log',errors='ignore'):
    p=l.replace(',',' ').split()
    if len(p)==3 and p[2].isdigit(): k=p[0]+':'+p[1]; mem[k]=max(mem.get(k,0),int(p[2]))
if mem: print('GPU mem peak MiB max %d min %d (n=%d)'%(max(mem.values()),min(mem.values()),len(mem)))
if M=='rt':
    def data(f):
        r=open(f,'rb').read(); i=r.find(b'<par_end>'); return r[i:] if i>=0 else r
    for h in sorted(glob.glob(D+'/a/*.hst')):
        b=os.path.join(D,'b',os.path.basename(h)); la=[l for l in open(h) if not l.startswith('#')]
        lb=[l for l in open(b) if not l.startswith('#')] if os.path.exists(b) else []
        print('rt',os.path.basename(h),'last row', 'same' if lb and la[-1]==lb[-1] else 'DIFF')
    last={}
    for f in sorted(glob.glob(D+'/a/bin/*.bin')): last[f.split('.')[-3]]=f
    for k,f in sorted(last.items()):
        g=os.path.join(D,'b','bin',os.path.basename(f)); print('rt',os.path.basename(f),'same' if os.path.exists(g) and data(f)==data(g) else 'DIFF/missing')
    for x in glob.glob(D+'/b/run*.log'): print('CW b',x[-8:],sum('CLEAN WARNING' in l for l in open(x,errors='ignore')))
