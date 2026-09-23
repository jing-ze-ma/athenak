#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/space2_3s/gate/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/space2_3s/gate/log.err.%j
#SBATCH -J m1space2
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# GPU cost of implicit_enthalpy = plm (runs_3s_space2): one binary, arms interleaved,
# repeat b in reverse order.  3-D He box 84x104x104, 4 blocks, 120 cycles,
# ms/cycle from cycle 20 to 120 (tools/summarize.py).
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/space2_3s
BIN=$W/new/build_gpu/src/athena
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
# implicit_enthalpy must be in the file (inp/box3d_plm = box3d_def + that one line)
declare -A A I
A[Eup]="";      I[Eup]=def
A[Eplm]="";     I[Eplm]=plm
A[Vup]="$V";    I[Vup]=def
A[Vplm]="$V";   I[Vplm]=plm
md5sum $BIN
arm () {
  local c=$1 r=$2; local name=T1_${c}_$r
  rm -rf $W/gate/runs/$name; mkdir -p $W/gate/runs/$name; cd $W/gate/runs/$name
  echo "#### $name $(hostname) $(date +%T)"
  srun -n 1 $BIN -i $W/gate/inp/box3d_${I[$c]}.athinput -d . -t 00:01:40 time/nlim=120 \
    ${A[$c]} > run.log 2> run.err
  echo "rc=$? $(date +%T)"; grep -E "Picard iterations|inner iterations|enthalpy" run.log | head -4
  cd $W/gate
}
ARMS="Eup Eplm Vup Vplm"
for c in $ARMS; do arm $c a; done
for c in $(echo $ARMS | tr ' ' '\n' | tac); do arm $c b; done
python3 $W/gate/tools/summarize.py T1_
