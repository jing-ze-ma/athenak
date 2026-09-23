#!/bin/bash
# is the Fx12r4 hybrid drift a Picard-count flip (sensitivity) or a bug?  the angle
# prototype (same ray split as G = 4) on the same input, and pure x1 splits
cd /viper/ptmp2/jinma/scscale_0923/cpu
O="time/tlim=3 output2/dt=1 output5/dt=1 output3/dt=1 problem/vpert=1.0e-3"
V="rad_m1/closure=vet_sc rad_m1/implicit_closure_lag=step"
F="rad_m1/vet_tensor=full"
X1_4="meshblock/nx1=21 rad_m1/implicit_partition=gather"
X12="meshblock/nx1=42 meshblock/nx2=16 rad_m1/implicit_partition=gather"
S=he_slab_mb_o3.athinput
DU="rad_m1/vet_dump=vet rad_m1/vet_dump_every=1"
./run.sh e_ref_Fx12r4 ref 4 $S $O $V $F $X12 $DU &
./run.sh e_ang_Fx12r4 new 4 $S $O $V $F $X12 $DU rad_m1/vet_mb_angles=true &
./run.sh e_hy4_Fx12r4 new 4 $S $O $V $F $X12 $DU rad_m1/vet_mb_agroup=4 &
./run.sh e_hy2_Fx12r4 new 4 $S $O $V $F $X12 $DU rad_m1/vet_mb_agroup=2 &
./run.sh e_ref_Ux1r4 ref 4 $S $O $V $X1_4 $DU &
./run.sh e_hy2_Ux1r4 new 4 $S $O $V $X1_4 $DU rad_m1/vet_mb_agroup=2 &
./run.sh e_hy4_Ux1r4 new 4 $S $O $V $X1_4 $DU rad_m1/vet_mb_agroup=4 &
./run.sh e_ang_Ux1r4 new 4 $S $O $V $X1_4 $DU rad_m1/vet_mb_angles=true &
wait
for d in e_*; do echo "$d: $(grep -E 'Picard iterations' $d/log.txt | cut -c1-140)"; done
python3 cmp_planes.py e_ref_Fx12r4 e_ang_Fx12r4 e_hy4_Fx12r4 e_hy2_Fx12r4 | uniq
python3 cmp_planes.py e_ref_Ux1r4 e_ang_Ux1r4 e_hy4_Ux1r4 e_hy2_Ux1r4 | uniq
