#!/bin/bash
# pagesync.sh <rsgE|hggl> [status]   (10-21) one-command refresh of a global-run results page, token-cheap.
#   no 2nd arg: submits 3-D movie job (incremental) -> page build job (afterany) -> stage into ~/ATHENAK/pages/<name>_glob/ ; prints 1 line.
#   status    : prints 1 line (job state, or READY + page time/sim time/sizes) -> then publish with the call in the memory file.
# RSG E: url LBASg2iK2ArcLY1V2yBCAR ; HG global: url LNtjtaL34oUhgcCRNp5V2M
J=/resnick/groups/carnegie_poc/jingze; N=$1
case $N in rsgE) PG=$J/rsgE_page; V=rsgE_glob; RD=$J/rsg_E_1020;; hggl) PG=$J/hgglob_page; V=hggl_glob; RD=$J/hgglob_1021;; *) echo "usage: pagesync.sh rsgE|hggl [status]"; exit 1;; esac
DST=/resnick/home/jingze/ATHENAK/pages/${N}_glob; LOG=$PG/out/logs; mkdir -p $DST/movies $LOG
if [ "$2" = status ]; then
  q=$(squeue -h -u $USER -n vis3d_$V,pagesync_$N -o "%j:%T" | tr '\n' ' ')
  [ -n "$q" ] && { echo "BUSY $q"; exit 0; }
  [ -f $DST/STAGED ] || { echo "NOT STAGED (run pagesync.sh $N)"; exit 0; }
  echo "READY $(cat $DST/STAGED)"; exit 0
fi
[ -n "$(squeue -h -u $USER -n vis3d_$V,pagesync_$N -o %i)" ] && { echo "already running"; exit 0; }
rm -f $DST/STAGED
mkdir -p $RD/logsF; for f in $RD/linkF.*.out $RD/linkF3.*.out $RD/linkS*.out; do [ -f "$f" ] && ln -sf ../$(basename $f) $RD/logsF/link.$(basename $f | cut -d. -f2).out; done
A=$(sbatch --parsable -A carnegie_poc -p expansion -c 32 --mem=240G -t 02:00:00 -J vis3d_$V -o $J/vis3d/$V/job.%j.out \
     --export=ALL,WORKERS=24,VIEWS="single layer" $J/vis3d/render_cs.sh $V) || exit 1
B=$(sbatch --parsable -d afterany:$A -A carnegie_poc -p expansion -c 8 --mem=32G -t 03:00:00 -J pagesync_$N -o $LOG/pagesync.%j.out \
     --wrap "bash $PG/update.sh && rsync -a --delete $PG/out/publish/ $DST/ && \
       T=\$(grep -oE 't = [0-9.e+]+ s = <b>[0-9.]+ d' $DST/index.html | sed 's/<b>//' | head -1); \
       echo \"\$(date '+%m-%d %H:%M') sim \$T | html \$(du -k $DST/index.html | cut -f1) kB, movies \$(du -ck $DST/movies/*.mp4 | tail -1 | cut -f1) kB\" > $DST/STAGED") || exit 1
echo "submitted 3-D $A -> page $B; check: bash $0 $N status"
