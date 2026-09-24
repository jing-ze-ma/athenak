#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/vetcol2_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/vetcol2_0924/gpu/log.err.%j
#SBATCH -J vetcol2_tune
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# team build tuning: vet_col_team_size x vet_col_chunk on the He wedge grid, same binary,
# interleaved, 2 repeats; Eddington and the one-thread kernel as references
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/vetcol2_0924; G=$W/gpu; D=/viper/ptmp2/jinma/defaults_0923/cpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <ref|new> <prob> <np> <args...>
  local t=$1 c=$2 p=$3 np=$4; shift 4
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_${p}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
H=$W/inp2/hewedge.athinput
A="rad_m1/closure=vet_col"
for r in 1 2; do
  run tn_edd_$r new none 2 -i $H rad_m1/closure=eddington
  run tn_c1_$r new none 2 -i $H $A rad_m1/vet_col_team=false
  for ts in 0 64 128; do for lc in 8 16 32; do
    run tn_t${ts}_l${lc}_$r new none 2 -i $H $A rad_m1/vet_col_team_size=$ts rad_m1/vet_col_chunk=$lc
  done; done
done
for t in $(ls $G | grep "^tn_"); do
  echo "== $t: $(grep -h 'cpu time used\|ms per build' $G/$t/new/log.txt | sed 's/.*= //;s/.*fenced), //' | tr '\n' ' ')"
done
echo TUNE DONE
