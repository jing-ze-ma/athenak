# one <bin> <pkg> <runroot> <tag> <np> <gpu 0|1|all>: one run inside the job, rocm-smi VRAM log
one() {
  local B=$1 K=$2 R=$3 tag=$4 np=$5 g=$6 d=$3/$4
  [ -d $d ] && mv $d $d.old.$(date +%m%d%H%M%S); mkdir -p $d; cd $d || return
  cp $K/decks/$tag.athinput in.athinput
  if [ $g = all ]; then
    cmd="srun --overlap -N1 -n $np --ntasks-per-node=$np --gres=gpu:2 -c 12 --cpu-bind=cores -J p16_$tag"
  else
    cmd="srun --overlap -N1 -n 1 --ntasks-per-node=1 --gres=gpu:2 -c 24 --cpu-bind=none -J p16_$tag env ROCR_VISIBLE_DEVICES=$g"
  fi
  echo "timeout 6m $cmd $B -i in.athinput -d ." > cmd.txt
  ( while true; do echo "$(date +%T) $(rocm-smi --showmeminfo vram --csv 2>/dev/null | grep -E '^card' | tr '\n' ' ')" >> mem.log; sleep 1; done ) &
  local mp=$!
  local t0=$(date +%s.%N)
  timeout 6m $cmd $B -i in.athinput -d . < /dev/null > run.log 2>&1; local rc=$?
  kill $mp 2>/dev/null
  echo "wall $(awk -v a=$t0 -v b=$(date +%s.%N) 'BEGIN{printf "%.2f", b-a}')" >> run.log; echo "rc=$rc" >> run.log
  echo "$(date +%T) $(basename $R) $tag gpu=$g rc=$rc fatal/nan=$(grep -c -E '### FATAL|(^|[^A-Za-z_])-?(nan|NaN)([^A-Za-z_]|$)|DIVERGED' run.log) warn=$(grep -c 'CLEAN WARNING' run.log) retry=$(grep -ci 'DT_RETRY' run.log)"
}
