#!/bin/bash
#SBATCH -A carnegie_poc -p expansion -N 1 -c 32 --mem=200G -t 01:30:00 -J pgrered -o /resnick/groups/carnegie_poc/jingze/rereduce.%j.out
# rereduce_job.sh rsgE|hggl K NS (10-24): re-reduce every dump to reduction VER 3 (gas+rad enthalpy flux, KE flux), shard K of NS, 32-wide
J=/resnick/groups/carnegie_poc/jingze; N=$1; K=$2; NS=$3
case $N in rsgE) PG=$J/rsgE_page;; hggl) PG=$J/hgglob_page;; esac
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
cd $PG/out
ls $BINDIR/*.hydro_w.*.bin | sed -E 's/.*\.([0-9]+)\.bin$/\1/' | awk -v k=$K -v n=$NS 'NR%n==k' \
  | xargs -r -P 32 -I{} env ONLY={} $PY $SRC/reduce.py $RUN data > /dev/null 2>&1
echo "$N shard $K/$NS done $(date +%T); ver3 files: $($PY -c "import glob,numpy as np;print(sum(int(np.load(f)['ver'])>=4 for f in glob.glob('data/red_*.npz')))")"
