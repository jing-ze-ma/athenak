#!/bin/bash
# compare R2 vs n48/R, RV vs CV (also R2 vs C, RV vs R for context): hst rows + bin DATA (header excluded)
P=/resnick/groups/carnegie_poc/jingze/rsg_g25_prod; S=$P/det/rdiff_1019; N=$P/det/n48
BC=/resnick/groups/carnegie_poc/jingze/hegiant_large_1014/bincmp.py
hc() { python3 - "$@" <<'PY'
import sys
def rows(p):
    d={}
    for l in open(p):
        if l.startswith('#'): continue
        d[l.split()[0]]=l
    return d
a,b=rows(sys.argv[1]),rows(sys.argv[2]); ks=sorted(set(a)&set(b),key=float)
bad=[k for k in ks if a[k]!=b[k]]
cyc=lambda k: round(float(k)/328.5714285714286)
print(sys.argv[3], 'IDENT' if not bad and ks else 'DIFF', len(ks),'rows', ('first diff t=%s (cycle %d)'%(bad[0],cyc(bad[0]))) if bad else '')
PY
}
pair() { A=$1; B=$2; L=$3; for h in rsg.hydro.hst rsg.user.hst; do [ -f $A/$h ] && [ -f $B/$h ] && hc $A/$h $B/$h "$L $h"; done
  for f in $A/bin/*.bin; do n=$(basename $f); [ -f $B/bin/$n ] && python3 $BC $f $B/bin/$n | sed "s/^/$L /"; done; }
pair $S/R2 $N/R "R2-R"
pair $S/RV $S/CV "RV-CV"
pair $S/R2 $N/C "R2-C"
pair $S/CV $N/C "CV-C"
