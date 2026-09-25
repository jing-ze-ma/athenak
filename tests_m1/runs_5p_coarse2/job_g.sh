#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/coarse2_0925/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/coarse2_0925/runs/log.err.%j
#SBATCH -J gcG
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-coarse2 GPU gates: (1) default bitwise base 038148c6 vs new (3 Cartesian cases of
# runs_5i_sphhalo + the wedge 1 GPU); (2) T-S4 vet_col atmosphere, 800 cycles: mg_gc vs
# the default, 1 GPU (4x4) and 2 GPUs (8x8), 1 x 1 and 2 x 2 bands.  Fresh run dirs.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/coarse2_0925; G=$W/runs/bw; D=/viper/ptmp2/jinma/defaults_0923/cpu
C=/viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_5i_sphhalo/scripts/cmpdir.py
TC=/viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_5i_sphhalo/scripts/tabcmp.py
NX=$W/bin/athena_g4_box_convection_gpu; RX=$W/bin/athena_base_box_convection_gpu
NW=$W/bin/athena_g4_none_gpu; RW=$W/bin/athena_base_none_gpu
md5sum $NX $RX $NW $RW
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() { local d=$G/$1 np=$2 b=$3; shift 3
  mkdir -p $d; cd $d; echo "== $1 $(date +%T)"
  timeout 240 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
    $b -d $d -t 00:03:00 "$@" > run.log 2> run.err
  echo "rc=$?" >> run.log; }
for c in ref new; do
  X=$NX; Y=$NW; [ $c = ref ] && X=$RX && Y=$RW
  run slab1/$c 1 $X -i $D/slab2d_old.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2/$c 2 $X -i $D/box3d_nd.athinput $B3
  run boxvet2/$c 2 $X -i $D/box3d_nd.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
  run wed1/$c 1 $Y -i /viper/ptmp2/jinma/mgf_0925/inp/hewedge_d.athinput time/nlim=20
done
for t in slab1 box2 boxvet2 wed1; do
  echo "== $t: $(python3 $C $G/$t/ref $G/$t/new) $(tail -n1 $G/$t/ref/run.log) $(tail -n1 $G/$t/new/run.log)"
done
R=$W/runs/ts4
IS=/viper/ptmp2/jinma/sphhalo_0924/inp/sph_atm_vc.athinput; IM=$W/inp/sph_atm_vc_m.athinput
Q="mesh/nx2=8 mesh/nx3=8"
GC="rad_m1/implicit_precond=mg_gc rad_m1/implicit_mg_levels=3 rad_m1/implicit_bcg_rho_direct=true"
B22="rad_m1/implicit_gc_bands2=2 rad_m1/implicit_gc_bands3=2"
runs() {  local n=$1 np=$2 i=$3; shift 3
  local d=$R/$n; mkdir -p $d; cd $d; echo "#### $n $(date +%T)"
  timeout 200 srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
      $NW -i $i -d . "$@" > run.log 2> run.err
  echo "rc=$? $(date +%T)"; }
runs def1 1 $IS
runs gc1 1 $IM $GC
runs gcb1 1 $IM $GC $B22
runs def2 2 $IS $Q
runs gc2 2 $IM $Q $GC
runs gcb2 2 $IM $Q $GC $B22
for p in "def1 gc1" "def1 gcb1" "def2 gc2" "def2 gcb2"; do set -- $p
  echo "T-S4 $2 vs $1: $(python3 $TC $R/$1 $R/$2) nonconv $(grep -o 'NON-CONVERGED=[0-9.e+-]*' $R/$2/run.log | tail -1)"
done
