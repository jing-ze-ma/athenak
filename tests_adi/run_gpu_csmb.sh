#!/bin/bash -l
#SBATCH -o gpu_csmb.out.%j
#SBATCH -e gpu_csmb.err.%j
#SBATCH -J adi_csmb
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# MORE THAN ONE MeshBlock PER PANEL under rad_ang_solver = adi, on one node.
# Each arm is run on 1 rank and on 2 ranks with the SAME MeshBlock layout, so what is
# compared is purely the rank decomposition: the blocks of a panel then live on different
# ranks and the open-chain interface gather runs over MPI.  The single-block-per-panel
# reference is run too, so the block decomposition can be separated from the rank one.
# The dumps are measured inside the job and deleted (inode quota).
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
B=$W/build_gpu_cs_test/src/athena
I=$W/inputs/tests/cubed_sphere_raddiff.athinput
# run <tag> <nranks> <lharm> <args...>
run(){ d=$1; nr=$2; lh=$3; shift 3
  rm -rf $d; mkdir -p $d
  (cd $d && srun -n $nr $B -i $I hydro/rad_ang_verbose=true "$@" > log.txt 2>&1)
  echo -n "$d  "
  python3 measure_cs.py $d $lh 2>/dev/null | sed "s/^[^ ]* *//"
  grep -o "|sum V de|/sum V|de| = [0-9.e+-]*" $d/log.txt | tail -1
  rm -rf $d/bin
}
for lh in 2 6; do
  # (A) n = 32 per panel edge: 1x1 and 2x2 MeshBlocks per panel, 1 and 2 ranks
  for s in sts adi; do
    for cfg in "32 1" "16 1" "16 2"; do
      set -- $cfg; mb=$1; nr=$2
      run C32_${s}_mb${mb}_r${nr}_l${lh} $nr $lh \
          mesh/nx2=32 mesh/nx3=32 meshblock/nx2=$mb meshblock/nx3=$mb \
          time/cfl_number=0.0075 time/nlim=80 problem/lharm=$lh \
          hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
    done
  done
done
# (B) n = 64 per panel edge with 4x4 MeshBlocks per panel (the production block count),
# 1 and 2 ranks, l = 2 only
for s in sts adi; do
  for cfg in "64 1" "16 1" "16 2"; do
    set -- $cfg; mb=$1; nr=$2
    run C64_${s}_mb${mb}_r${nr}_l2 $nr 2 \
        mesh/nx2=64 mesh/nx3=64 meshblock/nx2=$mb meshblock/nx3=$mb \
        time/cfl_number=0.0075 time/nlim=80 \
        hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
  done
done
echo GPU_CSMB_DONE
