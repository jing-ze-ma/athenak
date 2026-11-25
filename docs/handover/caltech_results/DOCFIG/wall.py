"""Wall time per star on one node of 4 H200, from log.N 'elapsed=' lines, chunks.txt and sacct step times. read only."""
import re, glob, os, subprocess, time, datetime as dt
J='/resnick/groups/carnegie_poc/jingze/'
R=dict(bsg=('prod_c1_1124/bsg','bsg3d',4.9248e6,[('log.0','7333035.52')]),m35mms=('prod_c1_1124/m35mms','m35mms',8.64e5,[('log.0','7333073.143'),('log.1','7333073.145')]),
 m35zams=('prod_c1_1124/m35zams','m35zams',8.64e5,[('log.0','7333073.144'),('log.1','7333073.146')]),t9l62=('prod_c4_1124/t9l62','t9l62',2.682e6,[('log.0','7333036.150'),('log.1','7333036.151')]),t9mlt=('prod_mlt_1125/t9mlt','t9l62',2.682e6,[('log.0','7333072.428')]))
CELLS=dict(bsg=256**3,m35mms=336*256*256,m35zams=384*128*128,t9l62=512*256*512,t9mlt=512*256*512)
now=time.time(); print('read',time.strftime('%Y-%m-%d %H:%M:%S'))
def step(s):
    o=subprocess.run(['sacct','-j',s,'-n','-P','--format=JobID,State,Start,End,ElapsedRaw'],capture_output=True,text=True).stdout.strip().splitlines()
    for l in o:
        p=l.split('|')
        if p[0]==s: return p
    return None
for s,(d,b,tl,segs) in R.items():
    print('====',s,'cells',CELLS[s],'tlim',tl); tot=0; cyc_total=0; state=[]; chunks=open(J+d+'/chunks.txt').read().splitlines()
    for lg,sid in segs:
        el=[]; 
        for l in open(J+d+'/'+lg,errors='replace'):
            if l.startswith('elapsed='):
                m=re.match(r'elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)',l); el.append((float(m.group(1)),int(m.group(2)),float(m.group(3))))
        st=step(sid); elapsed_step=float(st[4]) if st and st[4] else None
        # s/cycle early/late: mean over the first and last 1000 cycles of the segment (excluding the first line)
        a=el[0]; z=el[-1]; n=len(el)
        def rate(i,j): return (el[j][0]-el[i][0])/(el[j][1]-el[i][1])
        e=[q for q in el if q[1]<=el[1][1]+1000]; early=rate(1,len(e)-1) if len(e)>3 else rate(0,n-1)
        l_=[k for k,q in enumerate(el) if q[1]>=z[1]-1000]; late=rate(l_[0],n-1)
        binm=[c for c in chunks if 'chunk %s '%re.sub('log.','',lg) in c and 'STAR' in c][0]; bn=re.search(r'athena_(\w+)',binm).group(0); md=re.search(r'md5=(\w{8})',binm).group(1)
        print('  %s: step %s state %s start %s end %s step elapsed %s s (%.3f h) | elapsed lines: first line t=%.1f s at cycle %d, last %.1f s at cycle %d time %.6e, rows %d | s/cycle early %.4f (c%d-%d) late %.4f (c%d-%d) whole %.4f | binary %s md5 %s'%(lg,sid,st[1] if st else '?',st[2] if st else '?',st[3] if st else '?',elapsed_step,(elapsed_step or 0)/3600,a[0],a[1],z[0],z[1],z[2],n,early,el[1][1],e[-1][1],late,el[l_[0]][1],z[1],rate(0,n-1),bn,md))
        tot+=elapsed_step or 0; state.append(st[1] if st else '?')
    H=J+d+'/'+b+'.user.hst'; last=open(H,'rb').read()[-4000:].decode('latin-1').splitlines()[-1].split(); t_end=float(last[0]); 
    rows=sum(1 for l in open(H,errors='replace') if not l.startswith('#') and l.strip()); print('  hst: last time %.8e, rows %d (cycles %d); sum of step elapsed over segments %.0f s = %.3f h; node-hours %.3f (1 node) ; model time per day of wall (tlim units): t %.4e s = %.3f d of model per day of wall: %.3f d/d'%(t_end,rows,rows-1,tot,tot/3600,tot/3600,t_end,t_end/86400,(t_end/86400)/(tot/86400)),'| complete' if t_end>=tl*(1-1e-9) else '| running: %.2f %% of tlim'%(100*t_end/tl), state)
