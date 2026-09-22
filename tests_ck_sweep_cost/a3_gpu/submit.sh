#!/bin/bash -l
#SBATCH -o log.out.%j
#SBATCH -e log.err.%j
#SBATCH -J cksw_a3
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# A2: five 300-cycle timing runs in ONE job, same restart/input as bench/prof_0922/plain.
#   r0  athena.base (pre-edit snapshot)       ck_spherical=true
#   r1  athena.opt  ck_sweep_cache=false      ck_spherical=true   (= r0 arithmetic)
#   r2  athena.opt  ck_sweep_cache=true       ck_spherical=true
#   r3  athena.opt  ck_sweep_cache=false      ck_spherical=false
#   r4  athena.opt  ck_sweep_cache=true       ck_spherical=false
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
B=/viper/u2/jinma/ATHENAK/bench/cksw_0922
N=4430131   # ncycle(4429831) + 300
run () {   # $1 tag  $2 binary  $3 sph  $4 cache
  d=$1; rm -rf $d; mkdir -p $d/bin $d/rst
  ( cd $d && srun -n 2 $2 \
      -r /viper/u2/jinma/ATHENAK/bench/cs_mhd_prod3/rst/dhj.00567.rst \
      -i ../deep_hot_jupiter.athinput \
      -t 00:02:30 \
      time/nlim=$N problem/ck_spherical=$3 problem/ck_sweep_cache=$4 \
      > run.log 2>run.err )
  echo "=== $1 sph=$3 cache=$4 : $(grep -h 'cpu time used' $d/run.log $d/run.err | tail -1)"
}
run r0 $B/athena.opt2 true  0
run r1 $B/athena.opt2 true  1
run r2 $B/athena.opt2 true  2
run r3 $B/athena.opt2 false 1
run r4 $B/athena.opt2 false 0
echo "---- summary ----"
grep -H "cpu time used\|zone-cycles" r*/run.log r*/run.err 2>/dev/null
