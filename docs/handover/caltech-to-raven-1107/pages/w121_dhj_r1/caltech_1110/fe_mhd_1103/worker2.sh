#!/bin/bash
# worker2.sh GPU|cpu [STOPFILE] : as worker.sh; NPROC from env (default 8); stops before a new claim when STOPFILE exists
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast; EOS=/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt
export NPROC=${NPROC:-8}
if [ "$1" = cpu ]; then X=gpu=0; export CUDA_VISIBLE_DEVICES=; else X=gpu=1; export CUDA_VISIBLE_DEVICES=$1; fi
STOP=$2
while true; do
  [ -n "$STOP" ] && [ -e $O/$STOP ] && break
  got=""
  while read pl d t; do mkdir $O/claim/$t 2>/dev/null && { got="$pl $d $t"; break; }; done < $O/queue.txt
  [ -n "$got" ] || break
  set -- $got; S=$(date +%s)
  nice -n 10 $FF/run_gpu.sh $1 $2 $3 $O/out $EOS $X
  echo "$3 $(( $(date +%s)-S )) s $(date +%T)"
  rm -f $O/out/spec_$3_pc?.npz $O/out_mx/spec_$3_pc?.npz
done; echo ALLDONE $(date +%T)
