#!/bin/bash -l
# restart gate: N cycles vs N/2 + restart + N/2, switch on (K=2, dev_halo=1)
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/launch_0923; exe=$W/bin/athena_new_boxcpu; G=$W/rstgate
SW="rad_m1/implicit_krylov_dev=2 rad_m1/implicit_krylov_dev_halo=1"
run() { # tag input N extra...
  t=$1; inp=$2; N=$3; shift 3; h=$((N/2))
  O="problem/vpert=1.0e-2 $SW"
  rm -rf $G/${t}_A $G/${t}_B; mkdir -p $G/${t}_A $G/${t}_B
  (cd $G/${t}_A && OMP_NUM_THREADS=1 $exe -i $G/$inp.athinput -d $G/${t}_A $O time/nlim=$N "$@" > log.txt 2>&1; echo rc=$? >> log.txt)
  (cd $G/${t}_B && OMP_NUM_THREADS=1 $exe -i $G/$inp.athinput -d $G/${t}_B $O time/nlim=$h "$@" > log1.txt 2>&1; echo rc=$? >> log1.txt
   r=$(ls rst/*.rst | tail -1); OMP_NUM_THREADS=1 $exe -r $r -d $G/${t}_B time/nlim=$N > log2.txt 2>&1; echo rc=$? >> log2.txt)
}
run bv box3d_def 30 mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 rad_m1/closure=vet_sc rad_m1/vet_tensor=full 2>/dev/null
run es slab2d_def 40
echo ALLDONE > $G/done.txt
