#!/bin/bash -l
# restart gate, implicit_halo_ovl_faces on, 2 ranks: N cycles vs N/2 + restart + N/2
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/sync_0924; exe=$W/bin/athena_new_boxcpu; G=$W/rstgate; N=30; h=15
O="problem/vpert=1.0e-2 mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 rad_m1/implicit_vimp_fold=true rad_m1/implicit_fast_kernels=true rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2 rad_m1/implicit_halo_mpi=true rad_m1/implicit_halo_overlap=true rad_m1/implicit_halo_ovl_faces=true"
M="mpirun -np 2 --oversubscribe --bind-to none"
rm -rf $G/A $G/B; mkdir -p $G/A $G/B
(cd $G/A && OMP_NUM_THREADS=1 $M $exe -i $W/cpu/box3d_be_x.athinput -d $G/A $O time/nlim=$N > log.txt 2>&1; echo rc=$? >> log.txt) &
(cd $G/B && OMP_NUM_THREADS=1 $M $exe -i $W/cpu/box3d_be_x.athinput -d $G/B $O time/nlim=$h > log1.txt 2>&1; echo rc=$? >> log1.txt
 r=$(ls rst/*.rst | tail -1); OMP_NUM_THREADS=1 $M $exe -r $r -d $G/B time/nlim=$N > log2.txt 2>&1; echo rc=$? >> log2.txt) &
wait
echo ALLDONE > $G/done.txt
