#!/bin/bash
# armc.sh TAG BIN RST NC KEYSFILE [EXTRA key=val ...]   (orion 2026-11-21: Caltech H200 port of our Raven/DeltaAI arm.sh)
# Runs one arm as an --overlap step inside a Caltech holder.
#   env (required): HOLDER (Slurm job id), NODE (one hostname), J (your root holding rsg_E_1020/ and hgglob_1021/), OUT (run root)
#   env (optional): NR (ranks, default 4 = 1 node x 4 H200), MODE (speed|rt), FINALRST (0|1, default 1), NH (rt: leg-b dump
#                   after NH cycles, default 1), KEEPBIN (1 = keep the .bin outputs), TL (srun time limit min, default 90),
#                   ENVX ("A=1,B=2" extra env for athena)
# speed: N0 -> N0+NC straight.  rt: a = N0 -> N0+NC straight; b = N0 -> N0+NH, write rst, restart from it -> N0+NC;
#        res.py compares the last hst rows and the final bins of a and b (bitwise expected).
TAG=$1; BIN=$2; RST=$3; NC=$4; KF=$5; shift 5; EX="$*"
: ${HOLDER:?} ${NODE:?} ${J:?} ${OUT:?}; ENVX=${ENVX//,/ }; NR=${NR:-4}; MODE=${MODE:-speed}
H=$(dirname "$(readlink -f "$0")")
case $(basename $RST) in rsg.01295.rst) N0=589788;; hegiant.00139.rst) N0=624961;; *) N0=${N0:?set N0};; esac
D=$OUT/$TAG; rm -rf $D; mkdir -p $D; cd $D
K=$(grep -v '^#' $KF | sed "s#@J@#$J#g" | tr '\n' ' ')
echo "TAG=$TAG BIN=$BIN md5=$(md5sum < $BIN | cut -c1-32) RST=$RST N0=$N0 NC=$NC MODE=$MODE HOLDER=$HOLDER NODE=$NODE NR=$NR ENVX=$ENVX EX=$EX" > info.txt
run() { local d=$1 n=$2; shift 2; mkdir -p $D/$d; cd $D/$d; t0=$(date +%s)
  OMPI_MCA_io=romio321 srun --overlap --jobid=$HOLDER -w $NODE -N 1 -n $NR --mpi=pmix -t ${TL:-90} -J c_$TAG \
    env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $ENVX $H/memwrap.sh $D $BIN "$@" $K time/nlim=$n $EX -d $D/$d \
    < /dev/null > $([ -f run.log ] && echo run2.log || echo run.log) 2>&1
  echo "rc=$? wall=$(( $(date +%s)-t0 ))" >> $D/$d/rc.txt; cd $D; }
E=$((N0+NC)); HH=$((N0+${NH:-1})); NOFIN=; [ "${FINALRST:-1}" = 0 ] && NOFIN=job/final_restart=false
run a $E -r $RST $NOFIN
if [ $MODE = rt ]; then
  run b $HH -r $RST; f=$(ls -t $D/b/rst/*.rst | head -1); echo "restart from $f" >> $D/b/rc.txt
  run b $E -r $f output5/dt=1e30 $NOFIN
fi
find $D -name "*.rst" -exec ls -l {} + > $D/rstsizes.txt 2>/dev/null; find $D -name "*.rst" -delete
python3 $H/res.py $D $MODE > $D/res.txt 2>&1; cat $D/*/rc.txt >> $D/res.txt
[ "${KEEPBIN:-0}" = 1 ] || find $D -name "*.bin" -delete
echo DONE >> $D/res.txt
