#!/bin/bash
# usage: mkrun.sh NAME TEMPLATE(L_Bd1|L_tc4|A_af1) "block/key=val ..."
# creates runs/NAME with athinput from the audit template + overrides
W=/raven/ptmp/jinma/athpp_fix_1016; A=/raven/ptmp/jinma/athpp_audit_1015/h2h/runs
name=$1; tpl=$2; shift 2
d=$W/runs/$name; mkdir -p $d
cp $A/$tpl/*.txt $d/ 2>/dev/null; rm -f $d/ledger*.txt $d/conv.txt
cp $A/$tpl/athinput $d/athinput
python3 - "$d/athinput" "$@" <<'PY'
import sys,re
f=sys.argv[1]; L=open(f).read().split('\n')
for ov in sys.argv[2:]:
    for kv in ov.split():
        blk,rest=kv.split('/'); key,val=rest.split('=')
        # find block
        bi=[i for i,l in enumerate(L) if l.strip()=='<%s>'%blk][0]
        ei=bi+1
        while ei<len(L) and not L[ei].startswith('<'): ei+=1
        done=False
        for i in range(bi+1,ei):
            if re.match(r'\s*%s\s*='%re.escape(key),L[i]): L[i]='%s = %s'%(key,val); done=True
        if not done: L.insert(ei if L[ei-1].strip() else ei-1,'%s = %s'%(key,val))
open(f,'w').write('\n'.join(L))
PY
echo $d
