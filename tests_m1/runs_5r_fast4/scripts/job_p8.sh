#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/fast4_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/fast4_0925/runs/log.err.%j
#SBATCH -J f4P8
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
W=/viper/ptmp2/jinma/fast4_0925; R=$W/runs/p8
for t in m4 m8; do md5sum $W/bin/athena_${t}_box_convection_gpu $W/bin/athena_${t}_none_gpu; done
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
IR=$W/inp/wedge_rst.athinput
KD="rad_m1/implicit_krylov_dev=2"; GN="rad_m1/implicit_gas_newton=true"
for r in 1 2 3; do
  for t in m4 m8; do
    BX=$W/bin/athena_${t}_box_convection_gpu; BW=$W/bin/athena_${t}_none_gpu
    run box1_${t}_$r 1 $BX $IB
    run wed1_${t}_$r 1 $BW $IR
    run box2_${t}_$r 2 $BX $IB $W2
    run wed2_${t}_$r 2 $BW $IR
  done
  BX=$W/bin/athena_m8_box_convection_gpu; BW=$W/bin/athena_m8_none_gpu
  run box1_m8kd_$r 1 $BX $IB $KD
  run box2_m8kd_$r 2 $BX $IB $W2 $KD
  run wed1_m8gn_$r 1 $BW $IW $GN
  run wed2_m8gn_$r 2 $BW $IW $GN
done
TM="rad_m1/implicit_timers=10"
for t in m4 m8; do
  BX=$W/bin/athena_${t}_box_convection_gpu; BW=$W/bin/athena_${t}_none_gpu
  run box1tm_$t 1 $BX $IB $TM
  run wed1tm_$t 1 $BW $IW $TM
  run box2tm_$t 2 $BX $IB $W2 $TM
  run wed2tm_$t 2 $BW $IW $TM
done
prof box1_m8 1 $W/bin/athena_m8_box_convection_gpu $IB
prof wed1_m8 1 $W/bin/athena_m8_none_gpu $IW
