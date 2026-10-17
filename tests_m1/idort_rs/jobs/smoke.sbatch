#!/bin/bash -l
#SBATCH -J rssmk
#SBATCH -p interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=120000
#SBATCH --time=00:40:00
#SBATCH -o /raven/ptmp/jinma/idortrs_1017/logs/smoke.%j.out
# idort-rs-1017: smoke (NC cycles) of every rsx deck + thread scaling on the biggest one
module purge; module load gcc/13 openmpi/5.0
W=/raven/ptmp/jinma/idortrs_1017; B=$W/bin/athena_cpu_b16b66db; NC=${NC:-30}; md5sum $B
S=$W/smoke/j.$SLURM_JOB_ID; mkdir -p $S
run() { r=$1; nt=$2; tag=$3; d=$S/$(echo $r | tr / _)_$tag; mkdir -p $d; cd $d
  t0=$(date +%s.%N); OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $B -i $W/runs/$r/in.athinput -d $d $(cat $W/runs/$r/args.txt) time/nlim=$NC > run.log 2>&1; rc=$?
  t1=$(date +%s.%N); echo "SMOKE $r nt=$nt rc=$rc $(echo "$t1 - $t0" | bc) s fatal=$(grep -c FATAL run.log) $(grep -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' run.log | tail -1) $(grep -o 'zone-cycles/cpu_second = .*' run.log)"; }
export -f run; export W B NC S
# all decks, 4 threads each, 16 at a time
awk '{print $1}' $W/runs/list_rsx.txt | xargs -P 16 -I{} bash -c 'run {} 4 a'
# scaling on M27 n2048 c0.2
for nt in 8 16 32; do run cart/rsx_m27/P_n2048_c0.2 $nt s$nt; done
echo SMOKE_DONE
