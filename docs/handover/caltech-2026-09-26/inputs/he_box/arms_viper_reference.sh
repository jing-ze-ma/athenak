# sourced by the job scripts: arm definitions for hebox_cfl_0926
W=/viper/ptmp2/jinma/hebox_cfl_0926
B=$W/bin/athena_gpu     # -> mr_0926 build of rt-integration 7d68dd70, md5 43929127
RST0=/viper/u2/jinma/ATHENAK/bench/m1_vet3dcfl_0923/C15/rst/m1slab.00008.rst   # t=3800, read in place
T0=3800
COMMON="-i $W/inp/box.athinput time/nlim=-1
 output5/dt=250 output5/last_time=$T0 output2/dt=1e9 output3/dt=2000 output3/last_time=$T0"
declare -A K
K[R3]="time/cfl_number=0.3"
K[Pn]="time/cfl_number=0.3"                       # run on 2 GPUs (round-off noise reference)
K[Pu]="time/cfl_number=0.30000000000000004"        # 1-ulp dt kick (1-GPU noise reference)
K[H9]="time/cfl_number=0.9"
K[H9S2]="time/cfl_number=0.9 rad_m1/vet_sc_every=2"
# run_arm <arm> <np> <gpu0> <tlim> <walllimit hh:mm:ss> <rundir>; restarts from the newest rst
# in <rundir>/rst if there is one, else from RST0
run_arm() { local a=$1 np=$2 g0=$3 tl=$4 wl=$5 d=$6
  mkdir -p $d; cd $d
  local prev; prev=$(ls -t $d/rst/*.rst 2>/dev/null | head -1)
  [ -n "$prev" ] || prev=$RST0
  echo "#### $a start from $prev tlim=$tl $(date +%T)"
  GPU0=$g0 srun --export=ALL -n $np --overlap --exact --cpus-per-task=24 bash -c \
    'export ROCR_VISIBLE_DEVICES=$((SLURM_LOCALID + GPU0)); exec "$0" "$@"' \
    $B -r $prev $COMMON ${K[$a]} -d $d -t $wl time/tlim=$tl >> $d/run.log 2>> $d/run.err
  echo "#### $a rc=$? $(date +%T)"; }
