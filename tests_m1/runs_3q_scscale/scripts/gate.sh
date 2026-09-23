#!/bin/bash
# CPU correctness gates for sc-scale: ref = HEAD 5cd0a86d, new = sc-scale
cd /viper/ptmp2/jinma/scscale_0923/cpu
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
DU="rad_m1/vet_dump=vet rad_m1/vet_dump_every=5"
FU=rad_m1/vet_mb_mom_fuse=true
NA=rad_m1/vet_mb_agg=false
H2=rad_m1/vet_mb_agroup=2; H4=rad_m1/vet_mb_agroup=4
if [ "$1" != cmp ]; then
for w in ref new; do
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
done
wait
# exact options (hst must be identical to ref)
./run.sh na_Ux2r4 new 4 $S $O $V $X2_4 $NA &
./run.sh na_T3r4 new 4 $T $T3 $NA &
./run.sh fu_Ux2b2r2 new 2 $S $O $V $X2_2 $FU &
./run.sh fu_Ux2r4 new 4 $S $O $V $X2_4 $FU &
./run.sh fu_Ux1r4 new 4 $S $O $V $X1_4 $FU &
./run.sh fu_Fx12r4 new 4 $S $O $V $F $X12 $FU &
./run.sh fu_Uoddr3 new 3 $S $O $V $XODD $FU &
./run.sh fu_T3r4 new 4 $T $T3 $FU &
./run.sh fu_T3r1 new 0 $T $T3 $FU &
./run.sh fu_T3oddr4 new 4 $T $T3 $T3O $FU &
./run.sh fuh3_T3r4 new 4 $T $T3 $FU rad_m1/vet_mb_halo=3 &
./run.sh fuh2o_T3r4 new 4 $T $T3 $FU rad_m1/vet_mb_halo=2 rad_m1/vet_mb_overlap=true &
./run.sh h3_T3r4 new 4 $T $T3 rad_m1/vet_mb_halo=3 &
./run.sh nah3o_T3r2 new 2 $T $T3 $NA rad_m1/vet_mb_halo=3 rad_m1/vet_mb_overlap=true &
wait
# hybrid (round-off): 2-D with plane dumps, 3-D hst
./run.sh dref_Ux2r4 ref 4 $S $O $V $X2_4 $DU &
./run.sh dref_Fx12r4 ref 4 $S $O $V $F $X12 $DU &
./run.sh hy2_Ux2r4 new 4 $S $O $V $X2_4 $DU $H2 &
./run.sh hy4_Ux2r4 new 4 $S $O $V $X2_4 $DU $H4 &
./run.sh hy2_Fx12r4 new 4 $S $O $V $F $X12 $DU $H2 &
./run.sh hy4_Fx12r4 new 4 $S $O $V $F $X12 $DU $H4 &
./run.sh hy2_T3r4 new 4 $T $T3 $H2 &
./run.sh hy4_T3r4 new 4 $T $T3 $H4 &
./run.sh hy2_T3r2 new 2 $T $T3 $H2 &
./run.sh hy2fuh3_T3r4 new 4 $T $T3 $H2 $FU rad_m1/vet_mb_halo=3 &
./run.sh hy2na_T3r4 new 4 $T $T3 $H2 $NA &
wait
fi
echo "== A. default arms, new vs ref (all files byte-identical)"
for a in U1 F1 Ux2r4 Fx2r4 Ux2r1 Ux2b2r2 Ux1r4 Fx12r4 Uoddr3 T3r4 T3r2 T3r1 T3oddr4; do
  ./compare.sh ref_$a new_$a; done
echo "== B. exact options vs ref (hst identical; bin/rst differ only in the echoed keys)"
for p in na_Ux2r4:Ux2r4 na_T3r4:T3r4 fu_Ux2b2r2:Ux2b2r2 fu_Ux2r4:Ux2r4 fu_Ux1r4:Ux1r4 \
         fu_Fx12r4:Fx12r4 fu_Uoddr3:Uoddr3 fu_T3r4:T3r4 fu_T3r1:T3r1 fu_T3oddr4:T3oddr4 \
         fuh3_T3r4:T3r4 fuh2o_T3r4:T3r4 h3_T3r4:T3r4 nah3o_T3r2:T3r2; do
  ./compare.sh ref_${p#*:} ${p%%:*}; done
echo "== C. hybrid vs ref (round-off)"
python3 cmp_planes.py dref_Ux2r4 hy2_Ux2r4 hy4_Ux2r4
python3 cmp_planes.py dref_Fx12r4 hy2_Fx12r4 hy4_Fx12r4
python3 hstdiff.py ref_T3r4 hy2_T3r4 hy4_T3r4 hy2_T3r2 hy2fuh3_T3r4 hy2na_T3r4
echo "== C2. hybrid variants against each other (exact options on top of G=2)"
./compare.sh hy2_T3r4 hy2fuh3_T3r4
./compare.sh hy2_T3r4 hy2na_T3r4
