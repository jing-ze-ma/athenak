#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/fast4_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/fast4_0925/runs/log.err.%j
#SBATCH -J f4S
#SBATCH --constraint="apu"
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-fast4 scaling at NP GPUs (argument; 2 per node): strong (box 84x208x208 = 16
# blocks, wedge 96x128x128 = 16 blocks) and weak (box 4 blocks of 52x52 per GPU, wedge
# 16 blocks of 32x32 per GPU), the m1-fast4 defaults (box input = vet_sc + mg 3, wedge = sp defaults) and the box with
# implicit_eos_cache_check_every = 10, same binary (TAG, default m8),
# interleaved, 2 reps.  Submit:
#   sbatch -p apudev -N 1 --ntasks-per-node=1 --gres=gpu:1 scripts/job_scale4.sh 1
#   sbatch -p apudev -N 1 scripts/job_scale4.sh 2
#   sbatch -p apu    -N 2 scripts/job_scale4.sh 4
#   sbatch -p apu    -N 4 scripts/job_scale4.sh 8
NP=$1
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/fast4_0925; TAG=${TAG:-m8}; R=$W/runs/s${NP}_$TAG
BX=$W/bin/athena_${TAG}_box_convection_gpu; BW=$W/bin/athena_${TAG}_none_gpu
md5sum $BX $BW
run() {  local n=$1 np=$2 b=$3 i=$4; shift 4
  local d=$R/$n; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 200 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
IB=$W/inp/box.athinput; IW=$W/inp/wedge.athinput
DX=1.0608776e7   # box dx2
BS="mesh/nx2=208 mesh/nx3=208 mesh/x2max=$(python3 -c "print(208*$DX)") mesh/x3max=$(python3 -c "print(208*$DX)")"
case $NP in 1) n2=104; n3=104;; 2) n2=208; n3=104;; 4) n2=208; n3=208;; 8) n2=416; n3=208;; esac
BWK="mesh/nx2=$n2 mesh/nx3=$n3 mesh/x2max=$(python3 -c "print($n2*$DX)") mesh/x3max=$(python3 -c "print($n3*$DX)")"
case $NP in 1) w2=128; w3=128; x3=1.5707963267948966;; 2) w2=128; w3=256; x3=3.141592653589793;;
  4) w2=128; w3=512; x3=6.283185307179586;; 8) w2=256; w3=512; x3=6.283185307179586;; esac
WWK="mesh/nx2=$w2 mesh/nx3=$w3 mesh/x3max=$x3"
EC="rad_m1/implicit_eos_cache_check_every=10"
for r in 1 2; do
  run boxS_k1_$r $NP $BX $IB $BS
  run boxS_ec_$r $NP $BX $IB $BS $EC
  run boxW_k1_$r $NP $BX $IB $BWK
  run boxW_ec_$r $NP $BX $IB $BWK $EC
  run wedS_k1_$r $NP $BW $IW
  run wedW_k1_$r $NP $BW $IW $WWK
done
# analysis: python3 scripts/tsum.py runs/s${NP}_$TAG 10 30
