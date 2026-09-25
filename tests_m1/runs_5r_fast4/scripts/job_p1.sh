#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/fast4_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/fast4_0925/runs/log.err.%j
#SBATCH -J f4P1
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-fast4 step 1: baseline profile at the current defaults (binary TAG, default t1)
# box = vet_sc + mg 3, wedge = vet_col + mg_gc 1; 1 and 2 GPUs; timing (2 reps),
# fenced code timers (implicit_timers=10), rocprofv3 kernel traces
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/fast4_0925; TAG=${TAG:-t1}; R=$W/runs/p1_$TAG
BX=$W/bin/athena_${TAG}_box_convection_gpu; BW=$W/bin/athena_${TAG}_none_gpu
md5sum $BX $BW
IB=$W/inp/box.athinput; IW=$W/inp/wedge.athinput
W2="mesh/nx2=208 mesh/x2max=2.20662552e9"
run() {  local n=$1 np=$2 b=$3 i=$4; shift 4
  local d=$R/$n; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 200 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 time/nlim=60 rad_m1/implicit_picard_log=61 "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
prof() { local n=$1 np=$2 b=$3 i=$4; shift 4
  local d=$R/prof/$n; mkdir -p $d; cd $d; echo "#### prof $n"
  timeout 250 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --hip-runtime-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 time/nlim=40 "$@" > run.log 2> run.err
  echo "rc=$?"; }
for r in 1 2; do
  run box1_$r 1 $BX $IB
  run wed1_$r 1 $BW $IW
  run box2_$r 2 $BX $IB $W2
  run wed2_$r 2 $BW $IW
done
TM="rad_m1/implicit_timers=10"
run box1tm_1 1 $BX $IB $TM
run wed1tm_1 1 $BW $IW $TM
run box2tm_1 2 $BX $IB $W2 $TM
run wed2tm_1 2 $BW $IW $TM
prof box1 1 $BX $IB
prof wed1 1 $BW $IW
prof box2 2 $BX $IB $W2
prof wed2 2 $BW $IW
