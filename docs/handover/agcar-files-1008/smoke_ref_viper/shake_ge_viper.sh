#!/bin/bash
#SBATCH -J agcshkge
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --constraint=apu
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
#SBATCH -o /viper/ptmp2/jinma/lbv_1008/agcar/geos/shake/shake_ge.%j.out
# AG Car GENERAL-EOS 3-D shake-down DRAFT (agcar/geos) (not submitted), 1 node x 2 MI300A, 4 blocks 480x64x64 (2 per rank).
# usage: sbatch [-p apu --time=HH:MM:SS] shake.sh <A|B> <rundir> [extra keys]
#   smoke:  sbatch shake_ge.sh A /viper/ptmp2/jinma/lbv_1008/agcar/geos/shake/smokeA time/nlim=10
#   WALL (env, default 00:12:30) = the code's -t limit; set it ~10 min below --time for apu runs.
export SBATCH_EXPORT=NONE
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
BIN=/viper/ptmp2/jinma/builds/bin/athena_he_gpu72_2fbc6aa1_eft   # branch he-ic-eint-from-t
BIN_MD5=4c66621bc838501fbd948f2fed83866a
C=$1; D=$2; shift 2
case "$C" in
  A) IN=/viper/ptmp2/jinma/lbv_1008/agcar/geos/shake/agcar_shakeA_ge.athinput
     IC=/viper/ptmp2/jinma/lbv_1008/agcar/geos/icA/ic_agcar_A_ge.txt ;;
  B) IN=/viper/ptmp2/jinma/lbv_1008/agcar/geos/shake/agcar_shakeB_ge.athinput
     IC=/viper/ptmp2/jinma/lbv_1008/agcar/geos/icB/ic_agcar_B_ge.txt ;;
  *) echo "case must be A or B"; exit 1 ;;
esac
[ -n "$D" ] || { echo "no run dir"; exit 1; }
mkdir -p "$D" || exit 1
cd "$D" || exit 1
M=$(md5sum $BIN | cut -d' ' -f1)
[ "$M" == "$BIN_MD5" ] || { echo "md5 $M != $BIN_MD5"; exit 1; }
md5sum $BIN $IN $IC /viper/ptmp2/jinma/lbv_1008/agcar/geos/eos/eos_dump.txt /viper/ptmp2/jinma/lbv_1008/agcar/tables/*_tops_gs98_x0.36_z0.02.txt
echo "job $SLURM_JOB_ID case $C keys: $*  start $(date +%T)"
srun -n $SLURM_NTASKS --cpu-bind=cores $BIN -i $IN -d "$D" -t ${WALL:-00:12:30} "$@" > run.log 2>&1
echo "rc=$? $(date +%T)"
