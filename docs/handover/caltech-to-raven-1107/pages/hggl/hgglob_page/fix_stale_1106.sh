#!/bin/bash
# USER 11-06: re-reduce HG dumps whose reduction predates the bin (rewound runs rewrote the same index). Backups in out/stale_1106.
J=/resnick/groups/carnegie_poc/jingze; PG=$J/hgglob_page
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
cd $PG/out; S=$PG/out/stale_1106
(
 exec 7>$PG/out/.catchup.lock; flock 7; exec 9>$PG/out/.update.lock; flock 9
 $PY - > $S/list.txt <<'P'
import os,glob
B=os.environ['BINDIR']
for f in sorted(glob.glob('data/red_*.npz')):
    i=f[-9:-4]; b=f'{B}/hegiant.hydro_w.{i}.bin'
    if os.path.exists(b) and os.path.getmtime(b)>os.path.getmtime(f): print(i)
P
 echo "stale list: $(wc -l < $S/list.txt): $(paste -sd' ' $S/list.txt)"
 for i in $(cat $S/list.txt); do mv data/red_$i.npz $S/; for m in eq inner phot map; do [ -f frames/$m/$i.png ] && { mkdir -p $S/$m; mv frames/$m/$i.png $S/$m/; }; done; done
 xargs -a $S/list.txt -r -P 4 -I{} env ONLY={} nice -n 19 $PY $SRC/reduce.py $RUN data > $S/reduce.log 2>&1
 $PY - <<'P'
import glob,numpy as np,os
fs=sorted(glob.glob('data/red_*.npz')); t=np.array([float(np.load(f)['t']) for f in fs])/86400
print('n',len(fs),'nonmono',[(fs[k][-9:-4],t[k],t[k+1]) for k in np.where(np.diff(t)<=0)[0]])
S='stale_1106'
for i in open(S+'/list.txt').read().split():
    print(i,'old t',round(float(np.load(f'{S}/red_{i}.npz')['t'])/86400,3),'new t',round(float(np.load(f'data/red_{i}.npz')['t'])/86400,3))
P
)
NSH=4 bash $J/catchup_job.sh hggl
echo "fix_stale done rc=$? $(date +%T)"
