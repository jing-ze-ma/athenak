#!/bin/bash
# armc2.sh TAG BIN RST NC KEYSFILE [EXTRA key=val ...]   (orion H200 batch 4, 2026-11-20): armc.sh for 1..3 nodes
#   As armc.sh, but NODES (comma list, e.g. hpc-sm-01-17,hpc-sm-02-03) replaces NODE, and NR ranks are spread evenly
#   (NR/#nodes per node; 1 GPU per rank, KOKKOS_MAP_DEVICE_ID_BY=mpi_rank). All nodes of one arm must be in the same holder.
#   env: HOLDER NODES J OUT (required); NR (default 4 x #nodes), MODE (speed|rt), FINALRST, NH, KEEPBIN, TL, ENVX as armc.sh;
#        NSYS=1: rank 0 under "nsys profile -t cuda,nvtx,mpi" (output a/nsys_r0.nsys-rep; other ranks plain)
TAG=$1; BIN=$2; RST=$3; NC=$4; KF=$5; shift 5; EX="$*"
: ${HOLDER:?} ${NODES:?} ${J:?} ${OUT:?}; ENVX=${ENVX//,/ }
NN=$(echo $NODES | tr ',' '\n' | grep -c .); NR=${NR:-$((4*NN))}; MODE=${MODE:-speed}; PPN=$((NR/NN))
[ $((PPN*NN)) = $NR ] || { echo "NR=$NR not divisible by $NN nodes"; exit 2; }
H=$(dirname "$(readlink -f "$0")")
case $(basename $RST) in rsg.01295.rst) N0=589788;; hegiant.00139.rst) N0=624961;; *) N0=${N0:?set N0};; esac
D=$OUT/$TAG; rm -rf $D; mkdir -p $D; cd $D
K=$(grep -v '^#' $KF | sed "s#@J@#$J#g" | tr '\n' ' ')
echo "TAG=$TAG BIN=$BIN md5=$(md5sum < $BIN | cut -c1-32) RST=$RST N0=$N0 NC=$NC MODE=$MODE HOLDER=$HOLDER NODES=$NODES NR=$NR PPN=$PPN ENVX=$ENVX EX=$EX UCX=$(env | grep ^UCX_ | tr '\n' ' ')" > info.txt
W=$H/memwrap.sh; [ "${NSYS:-0}" = 1 ] && W=$H/nsyswrap.sh
run() { local d=$1 n=$2; shift 2; mkdir -p $D/$d; cd $D/$d; t0=$(date +%s)
  OMPI_MCA_io=romio321 srun --overlap --jobid=$HOLDER -w $NODES -N $NN -n $NR --ntasks-per-node=$PPN --mpi=pmix -t ${TL:-90} -J c_$TAG \
    env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $ENVX $W $D $BIN "$@" $K time/nlim=$n $EX -d $D/$d \
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
for h in $D/a/*.hst; do echo "hstmd5 $(basename $h) $(grep -v '^#' $h | md5sum | cut -c1-32) rows $(grep -vc '^#' $h)"; done >> $D/res.txt
[ "${KEEPBIN:-0}" = 1 ] || find $D -name "*.bin" -delete
echo DONE >> $D/res.txt
