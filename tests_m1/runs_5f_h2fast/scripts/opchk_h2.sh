#!/bin/bash -l
# implicit_op_check under hesdirk2 with the v4 defaults: cycles 1-3 (BE step + 2 x 2 stage solves)
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
W=/viper/ptmp2/jinma/h2fast_0924; E=$W/bin/athena_v4_box; G=$W/ginp; R=$W/cpu/opchk
BOX="box3d_h2d mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16"
VET="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
job() { local n=$1 np=$2; shift 2; local d=$R/$n; rm -rf $d; mkdir -p $d; set -- "$@"; local inp=$1; shift
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $E -i $G/$inp.athinput -d $d \
    problem/vpert=1.0e-2 time/nlim=3 rad_m1/implicit_op_check=-5 "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
job box1 1 $BOX
job box2h 2 $BOX rad_m1/implicit_halo_mpi=true
job boxv1 1 $BOX $VET
job slab1 1 slab2d_plm_vimp_h2d meshblock/nx2=32
job slab2h 2 slab2d_plm_vimp_h2d meshblock/nx2=16 rad_m1/implicit_halo_mpi=true
job slabv2 2 slab2d_plm_vimp_h2d meshblock/nx2=16 $VET
wait
for n in box1 box2h boxv1 slab1 slab2h slabv2; do L=$R/$n/log.txt
  echo "$n pass $(grep -c '^M1OPCHK solve [0-9]*: PASS' $L) fail $(grep -c '^M1OPCHK solve [0-9]*: FAIL' $L) $(grep -o 'rc=[0-9]*' $L | tail -1) worst $(grep '^  var' $L | awk '{print $4}' | sort -g | tail -1) fold=$(grep -o 'implicit_vimp_fold *= *[a-z]*' $L | head -1)"
done
echo OPCHK_DONE
