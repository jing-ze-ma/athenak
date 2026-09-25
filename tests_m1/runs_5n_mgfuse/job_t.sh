#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/mgf_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/mgf_0925/runs/log.err.%j
#SBATCH -J mgfT
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-mgfuse timing: rbgs_fwd / mg 3 unfused / mg 3 fused, same binary, interleaved,
# box vet_sc and wedge vet_col, 1 and 2 GPUs, + lin_tol 1e-8 arms (ms per iteration),
# + rocprofv3 kernel traces of the fused arms
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/mgf_0925; R=$W/runs/t
BX=$W/bin/athena_v1_box_convection_gpu; BW=$W/bin/athena_v1_none_gpu
md5sum $BX $BW
W2="mesh/nx2=208 mesh/x2max=2.20662552e9"
run() {  local n=$1 np=$2 b=$3 i=$4; shift 4
  local d=$R/$n; rm -rf $d; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 150 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -i $i -d . -t 00:02:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
IB=$W/inp/box_vsc_m.athinput; IW=$W/inp/hewedge_d.athinput
MU="rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=3 rad_m1/implicit_mg_fuse=false"
MF="rad_m1/implicit_precond=mg rad_m1/implicit_mg_levels=3 rad_m1/implicit_mg_fuse=true"
T8="rad_m1/implicit_lin_tol=1.0e-8"
for r in 1 2; do
  for a in def mu mf; do
    X=""; [ $a = mu ] && X=$MU; [ $a = mf ] && X=$MF
    run box1_${a}_$r 1 $BX $IB $X
    run box1_${a}T8_$r 1 $BX $IB $X $T8
    run box2_${a}_$r 2 $BX $IB $W2 $X
    run wed1_${a}_$r 1 $BW $IW $X
    run wed1_${a}T8_$r 1 $BW $IW $X $T8
    run wed2_${a}_$r 2 $BW $IW $X
  done
done
prof() { local n=$1 b=$2 i=$3; shift 3
  local d=$W/runs/prof_v1/$n; rm -rf $d; mkdir -p $d; cd $d; echo "#### prof $n"
  timeout 200 srun -n 1 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -i $i -d . -t 00:02:00 time/nlim=40 rad_m1/implicit_picard_log=41 "$@" > run.log 2> run.err
  echo "rc=$?"; }
prof box1_mu $BX $IB $MU
prof box1_mf $BX $IB $MF
prof wed1_mu $BW $IW $MU
prof wed1_mf $BW $IW $MF
