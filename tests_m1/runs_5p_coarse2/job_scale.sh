#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/coarse2_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/coarse2_0925/runs/log.err.%j
#SBATCH -J gcS
#SBATCH --constraint="apu"
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-coarse2 scaling at NP GPUs (argument; 2 per node): strong (box 84x208x208 = 16
# blocks, wedge 96x128x128 = 16 blocks) and weak (box 4 blocks of 52x52 per GPU, wedge
# 16 blocks of 32x32 per GPU), rbgs_fwd / mg 3 / mg_gc 1 / mg_gc 3, same binary,
# interleaved, 2 reps.  Submit:
#   sbatch -p apudev -N 1 --ntasks-per-node=1 --gres=gpu:1 scripts/job_scale.sh 1
#   sbatch -p apudev -N 1 scripts/job_scale.sh 2
#   sbatch -p apu    -N 2 scripts/job_scale.sh 4
#   sbatch -p apu    -N 4 scripts/job_scale.sh 8
NP=$1
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/coarse2_0925; R=$W/runs/s$NP
BX=$W/bin/athena_g4_box_convection_gpu; BW=$W/bin/athena_g4_none_gpu
md5sum $BX $BW
run() {  local n=$1 np=$2 b=$3 i=$4; shift 4
  local d=$R/$n; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 200 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
IB=$W/inp/box_vsc_m.athinput; IW=$W/inp/hewedge_d.athinput
DX=1.0608776e7   # box dx2
BS="mesh/nx2=208 mesh/nx3=208 mesh/x2max=$(python3 -c "print(208*$DX)") mesh/x3max=$(python3 -c "print(208*$DX)")"
case $NP in 1) n2=104; n3=104;; 2) n2=208; n3=104;; 4) n2=208; n3=208;; 8) n2=416; n3=208;; esac
BWK="mesh/nx2=$n2 mesh/nx3=$n3 mesh/x2max=$(python3 -c "print($n2*$DX)") mesh/x3max=$(python3 -c "print($n3*$DX)")"
case $NP in 1) w2=128; w3=128; x3=1.5707963267948966;; 2) w2=128; w3=256; x3=3.141592653589793;;
  4) w2=128; w3=512; x3=6.283185307179586;; 8) w2=256; w3=512; x3=6.283185307179586;; esac
WWK="mesh/nx2=$w2 mesh/nx3=$w3 mesh/x3max=$x3"
MG3="rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=3"
GC3="rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=3"
GC1="rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=1"
for r in 1 2; do
  for a in rb mg3 gc1 gc3; do
    X=""; [ $a = mg3 ] && X=$MG3; [ $a = gc3 ] && X=$GC3; [ $a = gc1 ] && X=$GC1
    run boxS_${a}_$r $NP $BX $IB $BS $X
    run boxW_${a}_$r $NP $BX $IB $BWK $X
    run wedS_${a}_$r $NP $BW $IW $X
    run wedW_${a}_$r $NP $BW $IW $WWK $X
  done
done
