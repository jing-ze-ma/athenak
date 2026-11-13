#!/usr/bin/env python3
"""blk.py deck [block=rad_m1] [args.txt]: print key=value of a block (deck then args overrides)."""
import sys,re
def read(deck, blk='rad_m1', args=None):
    d={}; b=None
    for ln in open(deck):
        s=ln.split('#')[0].strip()
        if not s: continue
        m=re.match(r'^<(\S+)>$',s)
        if m: b=m.group(1); continue
        if b==blk and '=' in s:
            k,v=s.split('=',1); d[k.strip()]=v.strip()
    if args:
        for a in open(args).read().split():
            if '=' in a and a.split('/',1)[0]==blk:
                k,v=a.split('/',1)[1].split('=',1); d[k]=v
    return d
if __name__=='__main__':
    blk=sys.argv[2] if len(sys.argv)>2 else 'rad_m1'
    for k,v in read(sys.argv[1],blk,sys.argv[3] if len(sys.argv)>3 else None).items(): print(f'{k}={v}')
