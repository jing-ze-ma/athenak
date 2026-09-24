#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/sphhalo_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/sphhalo_0924/gpu/log.err.%j
#SBATCH -J sphhalo
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-sphhalo GPU: Cartesian bitwise ref (d9583121) vs new, and the sp He wedge (vet_col,
# 2 GPUs) halo_mpi off vs on (default), same binary, interleaved, 3 repeats + profiles.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sphhalo_0924; G=$W/gpu; D=/viper/ptmp2/jinma/defaults_0923/cpu
md5sum $W/bin/*_gpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <dir> <np> <prof 0|1> <exe> <args...>
  local d=$G/$1 np=$2 p=$3 b=$4; shift 4
  rm -rf $d; mkdir -p $d; cd $d
  echo "== $1 $(date +%T)"
  if [ "$p" = 1 ]; then
    timeout 240 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --hip-runtime-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -d $d -t 00:03:00 "$@" > run.log 2> run.err
  else
    timeout 240 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -d $d -t 00:03:00 "$@" > run.log 2> run.err
  fi
  echo "rc=$?" >> run.log
}
NW=$W/bin/athena_new_none_gpu
# wedge timing, interleaved
for r in 1 2 3; do
  run w_off_$r 2 0 $NW -i $W/inp/hewedge_off.athinput time/nlim=40
  run w_on_$r 2 0 $NW -i $W/inp/hewedge.athinput time/nlim=40
done
run p_off 2 1 $NW -i $W/inp/hewedge_off.athinput time/nlim=40
run p_on 2 1 $NW -i $W/inp/hewedge.athinput time/nlim=40
# Cartesian bitwise
for c in ref new; do
  X=$W/bin/athena_${c}_box_convection_gpu
  run slab1/$c 1 0 $X -i $D/slab2d_old.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2/$c 2 0 $X -i $D/box3d_nd.athinput $B3
  run boxvet2/$c 2 0 $X -i $D/box3d_nd.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
done
for t in slab1 box2 boxvet2; do
  echo "== $t: $(python3 $W/scripts/cmpdir.py $G/$t/ref $G/$t/new) $(tail -n1 $G/$t/ref/run.log) $(tail -n1 $G/$t/new/run.log)"
done
for p in p_off p_on; do python3 $W/scripts/prof.py $G/$p/prof 0 40 10 > $G/$p.r0.sum 2>&1; done
echo GPU_DONE
