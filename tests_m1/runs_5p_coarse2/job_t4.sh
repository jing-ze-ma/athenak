#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/coarse2_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/coarse2_0925/runs/log.err.%j
#SBATCH -J gcP
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-coarse2 final timing (binary g4): rbgs_fwd / mg 3 / mg_gc 1, 2, 3, box and wedge,
# 1 and 2 GPUs, same binary, interleaved, 2 reps; lin_tol 1e-8 arms (ms/iteration);
# rocprofv3 traces of box mg_gc 3 and wedge mg_gc 1
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/coarse2_0925; R=$W/runs/t4
BX=$W/bin/athena_g4_box_convection_gpu; BW=$W/bin/athena_g4_none_gpu
md5sum $BX $BW
run() {  local n=$1 np=$2 b=$3 i=$4; shift 4
  local d=$R/$n; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 150 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -i $i -d . -t 00:02:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
IB=$W/inp/box_vsc_m.athinput; IW=$W/inp/hewedge_d.athinput
MG3="rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=3"
GC1="rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=1"
GC2="rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=2"
GC3="rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=3"
T8="rad_m1/implicit_lin_tol=1.0e-8"
W2="mesh/nx2=208 mesh/x2max=2.20662552e9"
for r in 1 2; do
  for a in rb mg3 gc1 gc2 gc3; do
    X=""; [ $a = mg3 ] && X=$MG3; [ $a = gc1 ] && X=$GC1; [ $a = gc2 ] && X=$GC2
    [ $a = gc3 ] && X=$GC3
    run box1_${a}_$r 1 $BX $IB $X
    run wed1_${a}_$r 1 $BW $IW $X
    run box2_${a}_$r 2 $BX $IB $W2 $X
    run wed2_${a}_$r 2 $BW $IW $X
    [ $r = 2 ] && run box1_${a}T8_$r 1 $BX $IB $X $T8
    [ $r = 2 ] && run wed1_${a}T8_$r 1 $BW $IW $X $T8
  done
done
prof() { local n=$1 b=$2 i=$3; shift 3
  local d=$W/runs/prof4/$n; mkdir -p $d; cd $d; echo "#### prof $n"
  timeout 200 srun -n 1 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -i $i -d . -t 00:02:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$?"; }
prof box1_gc3 $BX $IB $GC3
prof wed1_gc1 $BW $IW $GC1
