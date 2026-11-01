#!/bin/bash -l
# w121hfresh_1101: hydro fresh-start WASP-121b tests. Env: ARMS="name:deck:nlim:tlimwall ..." (sequential in one job)
#SBATCH -J w121hf1x
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --constraint=apu
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
#SBATCH -o /viper/ptmp2/jinma/w121hfresh_1101/logs/job.%j.out
#SBATCH -e /viper/ptmp2/jinma/w121hfresh_1101/logs/job.%j.err
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
D=/viper/ptmp2/jinma/w121hfresh_1101
BIN=$D/bin/athena_dhj_gpu72_18247755_hfresh
PKG=/viper/ptmp2/jinma/w121hfresh_1101/pkg_1x
cd $D || exit 1
md5sum $BIN
for a in $ARMS; do
  IFS=: read NAME DECK NL TW <<< "$a"
  R=$D/runs/$NAME
  mkdir -p $R && cd $R || exit 1
  ln -sfn $PKG/ck ck; ln -sfn $PKG/ic ic
  md5sum $D/deck/$DECK
  ( while true; do echo "T $(date +%s)"; grep -E "MemTotal|MemAvailable" /proc/meminfo; rocm-smi --showmeminfo vram gtt 2>/dev/null | grep -i "used"; sleep 15; done ) > $R/mem.log 2>&1 &
  MP=$!
  echo "$(date +%T) start $NAME"
  if [ "$NL" = "parse" ]; then
    $BIN -n -i $D/deck/$DECK > $R/parse.log 2>&1; echo "rc=$?"
  else
    srun -n 2 --cpu-bind=cores $BIN -i $D/deck/$DECK time/nlim=$NL -t $TW > $R/run.log 2>&1
    echo "$(date +%T) $NAME rc=$?"
  fi
  kill $MP
  cd $D || exit 1
done
