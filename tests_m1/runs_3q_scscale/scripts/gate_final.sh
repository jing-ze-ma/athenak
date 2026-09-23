#!/bin/bash
# final CPU gates with the final binary (bin/athena_new_cpu md5 0f9009c5): part A (new_* on the
# original inputs; ref_* from gate.sh), part B/C (gate2.sh), the Fx12 hybrid (gate3.sh)
cd /viper/ptmp2/jinma/scscale_0923/cpu
md5sum ../bin/athena_ref_cpu ../bin/athena_new_cpu
O="time/tlim=3 output2/dt=1 output5/dt=1 output3/dt=1 problem/vpert=1.0e-3"
V="rad_m1/closure=vet_sc rad_m1/implicit_closure_lag=step"
F="rad_m1/vet_tensor=full"
X2_2="meshblock/nx2=16"; X2_4="meshblock/nx2=8"
X1_4="meshblock/nx1=21 rad_m1/implicit_partition=gather"
X12="meshblock/nx1=42 meshblock/nx2=16 rad_m1/implicit_partition=gather"
XODD="mesh/nx1=81 meshblock/nx1=27 meshblock/nx2=16 rad_m1/implicit_partition=gather"
T3="mesh/nx2=32 mesh/x2max=3.39480832e8 mesh/nx3=32 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=6 output2/dt=0.5 output3/dt=0.5 rad_m1/closure=vet_sc"
T3O="mesh/nx1=83 meshblock/nx1=83"
S=he_slab_mb.athinput; T=he3d.athinput
w=new
./run.sh ${w}_U1 $w 0 $S $O $V &
./run.sh ${w}_F1 $w 0 $S $O $V $F &
./run.sh ${w}_Ux2r4 $w 4 $S $O $V $X2_4 &
./run.sh ${w}_Fx2r4 $w 4 $S $O $V $F $X2_4 &
./run.sh ${w}_Ux2r1 $w 0 $S $O $V $X2_4 &
./run.sh ${w}_Ux2b2r2 $w 2 $S $O $V $X2_2 &
./run.sh ${w}_Ux1r4 $w 4 $S $O $V $X1_4 &
./run.sh ${w}_Fx12r4 $w 4 $S $O $V $F $X12 &
./run.sh ${w}_Uoddr3 $w 3 $S $O $V $XODD &
./run.sh ${w}_T3r4 $w 4 $T $T3 &
./run.sh ${w}_T3r2 $w 2 $T $T3 &
./run.sh ${w}_T3r1 $w 0 $T $T3 &
./run.sh ${w}_T3oddr4 $w 4 $T $T3 $T3O &
wait
echo "== A. default arms, new vs ref (all files byte-identical)"
for a in U1 F1 Ux2r4 Fx2r4 Ux2r1 Ux2b2r2 Ux1r4 Fx12r4 Uoddr3 T3r4 T3r2 T3r1 T3oddr4; do
  ./compare.sh ref_$a new_$a; done
bash gate2.sh
bash gate3.sh
