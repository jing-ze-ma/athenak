#!/bin/bash -l
#SBATCH -J m1ovl_gate
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
#SBATCH -o /viper/ptmp2/jinma/m1int_0924/ovl_gate.%j.out
# Step 1 GPU gate: implicit_halo_overlap on 2 GPUs (30 cycles)
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/m1int_0924
md5sum $W/bin/athena_new_boxgpu $W/bin/athena_ovl_boxgpu
HM="rad_m1/implicit_halo_mpi=true"; PP="rad_m1/implicit_krylov_pipe=true"; OV="rad_m1/implicit_halo_overlap=true"
g() {  # <name> <new|ovl> <input> args...
  local n=$1 x=$W/bin/athena_$2_boxgpu inp=$3; shift 3
  local d=$W/gpu/$n; rm -rf $d; mkdir -p $d; cd $d
  srun -n 2 --ntasks-per-node=2 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
    $x -i $W/gpu/$inp -d $d meshblock/nx3=26 time/nlim=30 "$@" > log.txt 2>&1
  echo "$n rc=$? $(grep -o 'NON-CONVERGED=[0-9.e+-]*' log.txt | tail -1)"
}
g o_nhm new he3d_fast.athinput $HM
g o_hm  ovl he3d_fast.athinput $HM
g o_hmo ovl he3d_fast.athinput $HM $OV
g o_ph  ovl he3d_fast.athinput $HM $PP
g o_pho ovl he3d_fast.athinput $HM $PP $OV
g o_bhm  ovl box3d_plm_vimp_h2.athinput $HM
g o_bhmo ovl box3d_plm_vimp_h2.athinput $HM $OV
g o_bpho ovl box3d_plm_vimp_h2.athinput $HM $PP $OV
cd $W/gpu
for p in "o_nhm o_hm" "o_hm o_hmo" "o_ph o_pho" "o_bhm o_bhmo" "o_bhm o_bpho"; do
  set -- $p; cmp -s $1/m1slab.hydro.hst $2/m1slab.hydro.hst && cmp -s $1/m1slab.user.hst $2/m1slab.user.hst && echo "$2 vs $1 hst BITWISE" || python3 $W/scripts/cmp.py $W/gpu $1 $2
done
