#!/bin/bash -l
#SBATCH -J rsrun
#SBATCH -p interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=120000
#SBATCH --time=01:59:00
#SBATCH -o /raven/ptmp/jinma/idortrs_1017/logs/run.%j.out
# idort-rs-1017: rsx family full runs, b16b66db plain CPU, production keys (args.txt = FIN5.1 method keys + ROL2U).
# LIST lines "<dir> <threads>"; finished dirs (rc.txt = 0) skipped; unfinished ones RESTART from their newest rst
# (wall limit WL leaves time for the rst dump).  Out: runs/<dir>/out/{run*.log,rc.txt,tab/,rst/}
module purge; module load gcc/13 openmpi/5.0
W=/raven/ptmp/jinma/idortrs_1017; B=$W/bin/athena_cpu_b16b66db; md5sum $B
WL=${WL:-01:50:00}
one() { r=$1; nt=$2; d=$W/runs/$r/out; mkdir -p $d; cd $d
  { [ -f rc.txt ] && [ "$(cat rc.txt)" = 0 ]; } && { echo "SKIP $r"; return; }
  n=$(ls run*.log 2>/dev/null | wc -l); f=$(ls rst/*.rst 2>/dev/null | sort | tail -1)
  t0=$(date +%s)
  if [ -n "$f" ]; then OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $B -r $f -d $d -t $WL > run$n.log 2>&1; rc=$?
  else
    grep -q '^<output9>' $W/runs/$r/in.athinput || printf '\n<output9>\nfile_type = rst\ndt = 1.0e30\n' >> $W/runs/$r/in.athinput
    OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $B -i $W/runs/$r/in.athinput -d $d -t $WL $(cat $W/runs/$r/args.txt) > run$n.log 2>&1; rc=$?; fi
  tl=$(grep -E '^ *tlim' $W/runs/$r/in.athinput | awk '{print $3}')
  tt=$(grep -o 'time=[0-9.e+-]*' run$n.log | tail -1 | cut -d= -f2)
  done_=$(grep -c 'Terminating on time limit' run$n.log)
  [ $rc = 0 ] && [ "$done_" -ge 1 ] && echo 0 > rc.txt
  echo "RUN $r nt=$nt rc=$rc tlimhit=$done_ t=$tt/$tl $(( $(date +%s) - t0 )) s fatal=$(grep -c FATAL run$n.log) $(grep -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' run$n.log | tail -1) dtretry=$(grep -c '^DT_RETRY cycle' run$n.log)"; }
export -f one; export W B WL
grep -v '^#' ${LIST:?} | xargs -P ${P:-18} -L 1 bash -c 'one "$@"' _
echo RUN_DONE
