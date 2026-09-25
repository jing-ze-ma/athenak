#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/fast4_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/fast4_0925/runs/log.err.%j
#SBATCH -J f4P3
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-fast4: MR (implicit_mr_every 2, 4) and eos_cache_check_every=10 timings, box + wedge,
# 1 and 2 GPUs, binary TAG (default m1), interleaved, 2 reps; kokkos-trace profiles (labels)
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/fast4_0925; TAG=${TAG:-m4}; R=$W/runs/p3_$TAG
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
  timeout 250 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --kokkos-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 time/nlim=40 "$@" > run.log 2> run.err
  echo "rc=$?"; }
K2="rad_m1/implicit_mr_every=2"; K4="rad_m1/implicit_mr_every=4"
EC="rad_m1/implicit_eos_cache_check_every=10"
V2="rad_m1/vet_sc_every=2"
for r in 1 2 3; do
  run box1_k1_$r 1 $BX $IB
  run box1_ec_$r 1 $BX $IB $EC
  run box1_k2_$r 1 $BX $IB $K2
  run box1_k2ec_$r 1 $BX $IB $K2 $EC
  run box1_k4ec_$r 1 $BX $IB $K4 $EC
  run box1_v2ec_$r 1 $BX $IB $V2 $EC
  run wed1_k1_$r 1 $BW $IW
  run wed1_k2_$r 1 $BW $IW $K2
  run wed1_k4_$r 1 $BW $IW $K4
  run box2_k1_$r 2 $BX $IB $W2
  run box2_ec_$r 2 $BX $IB $W2 $EC
  run box2_k2ec_$r 2 $BX $IB $W2 $K2 $EC
  run box2_k4ec_$r 2 $BX $IB $W2 $K4 $EC
  run wed2_k1_$r 2 $BW $IW
  run wed2_k2_$r 2 $BW $IW $K2
  run wed2_k4_$r 2 $BW $IW $K4
done
TM="rad_m1/implicit_timers=10"
run box1tm_k1 1 $BX $IB $TM
run box1tm_k2ec 1 $BX $IB $TM $K2 $EC
run wed1tm_k1 1 $BW $IW $TM
run wed1tm_k2 1 $BW $IW $TM $K2
prof wed1 1 $BW $IW
