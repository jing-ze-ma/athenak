#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/mgf_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/mgf_0925/runs/log.err.%j
#SBATCH -J mgfprof
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
# usage: sbatch job_prof.sh <tag> <box binary> <wedge binary> [extra mg keys]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/mgf_0925; T=$1; BX=$2; BW=$3; shift 3
md5sum $BX $BW
run() { local n=$1 b=$2 i=$3; shift 3
  local d=$W/runs/prof_$T/$n; rm -rf $d; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 200 srun -n 1 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -i $i -d . -t 00:02:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
MG="rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=3"
run box1_mg3 $BX $W/inp/box_vsc_m.athinput $MG "$@"
run wed1_mg3 $BW $W/inp/hewedge_d.athinput $MG "$@"
run box1_def $BX $W/inp/box_vsc_m.athinput
