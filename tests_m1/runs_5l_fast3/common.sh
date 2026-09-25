# sourced by the job scripts
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/fast3
IB=$W/inp/box_vsc.athinput; IW=$W/inp/hewedge.athinput
W2="mesh/nx2=208 mesh/x2max=2.20662552e9"
SPL="rad_m1/implicit_op_split_red=true"
KD="rad_m1/implicit_krylov_dev=2 rad_m1/implicit_krylov_dev_halo=1"
PIPE="rad_m1/implicit_krylov_pipe=true"
run() {  # run <tag> <arm> <np> <prof 0|1> <bin tag> <prob> <input> [ov...]
  local t=$1 a=$2 np=$3 p=$4 b=$W/bin/athena_$5_$6 i=$7; shift 7
  local d=$W/runs/$t/$a; rm -rf $d; mkdir -p $d; cd $d
  echo "#### $t/$a $(date +%T)"
  if [ "$p" = 1 ]; then
    timeout 240 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --hip-runtime-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 "$@" > run.log 2> run.err
  else
    timeout 200 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $b -i $i -d . -t 00:03:00 "$@" > run.log 2> run.err
  fi
  echo "rc=$? $(date +%T)"
}
