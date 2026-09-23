#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/h2_3x/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/h2_3x/gpu/log.err.%j
#SBATCH -J h2g7b
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# runs_3x_hesdirk2 G0 (GPU, be bitwise vs 2c3c4178) + G7 (GPU cost hesdirk2 vs be): one
# binary, arms interleaved, repeat b reversed.  3-D He box 84x104x104, 4 blocks, plm+vimp,
# 120 cycles, ms/cycle from cycle 20 to 120.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/h2_3x
BIN=$W/new/build_gpu/src/athena
md5sum $BIN $W/base/build_gpu/src/athena
# --- G0 on the GPU
OV="time/nlim=40 output2/dt=2.0 output5/dt=2.0 output3/dt=1.0e9"
for b in base new; do
  d=$W/gpu/g0c_$b; rm -rf $d; mkdir -p $d; cd $d
  srun -n 1 $W/$b/build_gpu/src/athena -i $W/inp/box3d_plm_vimp.athinput -d . $OV > run.log 2>&1
  echo "g0 $b rc=$?"
done
cd $W/gpu
for f in $(cd g0c_new && find . -name "*.bin" -o -name "*.hst" | sort); do
  cmp -s g0c_new/$f g0c_base/$f && echo "SAME $f" || echo "DIFF $f"
done
# --- G7
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
declare -A A I
A[Ebe]="";  I[Ebe]=box3d_plm_vimp
A[Eh2]="";  I[Eh2]=box3d_plm_vimp_h2
A[Vbe]="$V"; I[Vbe]=box3d_plm_vimp
A[Vh2]="$V"; I[Vh2]=box3d_plm_vimp_h2
arm () {
  local c=$1 r=$2; local name=T8_${c}_$r
  rm -rf $W/gpu/runs/$name; mkdir -p $W/gpu/runs/$name; cd $W/gpu/runs/$name
  echo "#### $name $(hostname) $(date +%T)"
  srun -n 1 $BIN -i $W/inp/${I[$c]}.athinput -d . -t 00:01:40 time/nlim=120 \
    ${A[$c]} > run.log 2> run.err
  echo "rc=$? $(date +%T)"; grep -E "Picard iterations|hesdirk2:" run.log | head -3
  cd $W/gpu
}
ARMS="Ebe Eh2 Vbe Vh2"
for c in $ARMS; do arm $c a; done
for c in $(echo $ARMS | tr ' ' '\n' | tac); do arm $c b; done
for c in $ARMS; do arm $c c; done
python3 $W/gpu/tools/summarize.py T8_
