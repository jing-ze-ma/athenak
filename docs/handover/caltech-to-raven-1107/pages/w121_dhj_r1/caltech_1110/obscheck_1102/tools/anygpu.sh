#!/bin/bash
# anygpu.sh cmd... : run on the first free agent GPU (02-03 1/3, 01-11 1/3) via srun --overlap + lock
O=/resnick/groups/carnegie_poc/jingze/obscheck_1102
for ng in "5544871 hpc-sm-02-03 3" "5544871 hpc-sm-02-03 1" "5274196 hpc-sm-01-11 1" "5274196 hpc-sm-01-11 3"; do
  read J N G <<< "$ng"
  [ -d /resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj/lock_${N}_${G} ] && continue
  srun --jobid=$J --overlap -n1 -c8 -w $N $O/tools/gpurun.sh $N $G "$@"; rc=$?
  [ $rc -eq 3 -o $rc -eq 4 ] && continue
  exit $rc
done
echo "NO FREE GPU"; exit 5
