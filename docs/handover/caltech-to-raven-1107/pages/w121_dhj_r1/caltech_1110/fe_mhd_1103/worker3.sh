#!/bin/bash
# worker3.sh GPU : claims queue3.txt lines "planet dump tag [outdir] [key=val..]" and runs run_fast2.sh (1 cut, prep=gpu)
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast; EOS=/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt
export CUDA_VISIBLE_DEVICES=$1
while true; do
  bash $O/addm2.sh >> $O/logs/addm2.log
  got=""
  while read ln; do t=$(echo $ln | awk '{print $3}'); mkdir $O/claim3/$t 2>/dev/null && { got="$ln"; break; }; done < $O/queue3.txt
  [ -n "$got" ] || break
  set -- $got; pl=$1; d=$2; t=$3; od=${4:-$O/out}; shift 4 2>/dev/null || shift $#; S=$(date +%s)
  E=(); K=(); for a in "$@"; do case $a in FEHIRES_*=*) E+=("$a");; *) K+=("$a");; esac; done
  env "${E[@]}" nice -n 10 $FF/run_fast2.sh $pl $d $t $od $EOS "${K[@]}"
  echo "$t $(( $(date +%s)-S )) s $(date +%T)"
  rm -f $od/spec_${t}_pc6.npz ${od}_mx/spec_${t}_pc6.npz
done; echo ALLDONE $(date +%T)
