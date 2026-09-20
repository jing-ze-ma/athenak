#!/bin/bash -l
#SBATCH -o gpu_csmb2n.out.%j
#SBATCH -e gpu_csmb2n.err.%j
#SBATCH -J adi_csmb2n
#SBATCH -p apu
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:20:00
# TWO NODES, 4 ranks: panel seams AND intra-panel MeshBlock boundaries now cross a NODE
# boundary, so the open-chain interface gather and the module's own halo exchange both
# run over the interconnect.  Same arms as run_gpu_csmb.sh, plus the He-star FeCZ box
# production configuration compared bitwise against the 2-rank single-node run.
# Swap `-p apudev` for `-p apu` and raise --time if apudev refuses two nodes.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
B=$W/build_gpu_cs_test/src/athena
I=$W/inputs/tests/cubed_sphere_raddiff.athinput
run(){ d=$1; nr=$2; lh=$3; shift 3
  rm -rf $d; mkdir -p $d
  (cd $d && srun -n $nr $B -i $I hydro/rad_ang_verbose=true "$@" > log.txt 2>&1)
  echo -n "$d  "
  python3 measure_cs.py $d $lh 2>/dev/null | sed "s/^[^ ]* *//"
  grep -o "|sum V de|/sum V|de| = [0-9.e+-]*" $d/log.txt | tail -1
  rm -rf $d/bin
}
for lh in 2 6; do
  for s in sts adi; do
    run N32_${s}_mb16_r4_l${lh} 4 $lh \
        mesh/nx2=32 mesh/nx3=32 meshblock/nx2=16 meshblock/nx3=16 \
        time/cfl_number=0.0075 time/nlim=80 problem/lharm=$lh \
        hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
  done
done
for s in sts adi; do
  run N64_${s}_mb16_r4_l2 4 2 \
      mesh/nx2=64 mesh/nx3=64 meshblock/nx2=16 meshblock/nx3=16 \
      time/cfl_number=0.0075 time/nlim=80 \
      hydro/rad_implicit_ang=true hydro/rad_ang_solver=$s
done
# the He-star FeCZ box (the Cartesian production configuration) on 4 ranks over 2 nodes,
# for the bitwise comparison against the 2-rank single-node run of run_gpu_hebox.sh
BB=$W/build_gpu_box_convection/src/athena
IB=$PWD/he_box_w8_smoke.athinput
for nr in 2 4; do
  d=heboxr_n$nr; rm -rf $d; mkdir -p $d
  (cd $d && srun -n $nr $BB -i $IB time/nlim=30 > log.txt 2>&1)
done
python3 cmpbin.py heboxr_n2/bin heboxr_n4/bin
rm -rf heboxr_n2/bin heboxr_n4/bin
echo GPU_CSMB2N_DONE
