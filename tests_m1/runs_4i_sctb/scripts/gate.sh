#!/bin/bash
# CPU gates for vet_mb_tblock (m1-sctb).  Waves of <= 12 ranks.
cd /viper/ptmp2/jinma/sctb_0923/cpu
md5sum ../bin/athena_ref_cpu ../bin/athena_new_cpu
O="time/tlim=3 output2/dt=1 output5/dt=1 output3/dt=1 problem/vpert=1.0e-3"
V="rad_m1/closure=vet_sc rad_m1/implicit_closure_lag=step"
F="rad_m1/vet_tensor=full"
X2_2="meshblock/nx2=16"; X2_4="meshblock/nx2=8"
X1_4="meshblock/nx1=21 rad_m1/implicit_partition=gather"
X12="meshblock/nx1=42 meshblock/nx2=16 rad_m1/implicit_partition=gather"
XODD="mesh/nx1=81 meshblock/nx1=27 meshblock/nx2=16 rad_m1/implicit_partition=gather"
T3="mesh/nx2=32 mesh/x2max=3.39480832e8 mesh/nx3=32 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=6 output2/dt=0.5 output3/dt=0.5 rad_m1/closure=vet_sc"
S=he_slab_mb.athinput; T=he3d.athinput
TB3="rad_m1/vet_mb_tblock=3"
declare -A A
A[U1]="0 $S $O $V"
A[Fx2r4]="4 $S $O $V $F $X2_4"
A[Ux2r1]="0 $S $O $V $X2_4"
A[Ux2b2r2]="2 $S $O $V $X2_2"
A[Ux1r4]="4 $S $O $V $X1_4"
A[Fx12r4]="4 $S $O $V $F $X12"
A[Uoddr3]="3 $S $O $V $XODD"
A[T3r4]="4 $T $T3"
A[T3r2]="2 $T $T3"
A[T3r1]="0 $T $T3"
go() { # name which arm extra...
  local n=$1 w=$2 a=$3; shift 3; ./run.sh $n $w ${A[$a]} "$@"; }
# wave 1-2: default arms ref / new
go ref_Fx2r4 ref Fx2r4 & go new_Fx2r4 new Fx2r4 & go ref_Ux1r4 ref Ux1r4 & wait
go new_Ux1r4 new Ux1r4 & go ref_Fx12r4 ref Fx12r4 & go new_Fx12r4 new Fx12r4 & wait
go ref_T3r4 ref T3r4 & go new_T3r4 new T3r4 & go ref_Uoddr3 ref Uoddr3 & go ref_U1 ref U1 & wait
go new_Uoddr3 new Uoddr3 & go new_U1 new U1 & go ref_T3r2 ref T3r2 & go new_T3r2 new T3r2 & go ref_Ux2b2r2 ref Ux2b2r2 & wait
go new_Ux2b2r2 new Ux2b2r2 & go ref_T3r1 ref T3r1 & go new_T3r1 new T3r1 & go ref_Ux2r1 ref Ux2r1 & go new_Ux2r1 new Ux2r1 & wait
# tblock arms: the inputs spell out the new keys (he3d_tb: + vet_mb_halo = 3;
# he_slab_mb_tbd / he3d_tbd: without vet_mb_halo, so the default halo rule applies)
gt() { # name arm input-kind extra...
  local n=$1 a=$2 k=$3; shift 3
  local x=${A[$a]}; x=${x/he3d.athinput/he3d_$k.athinput}; x=${x/he_slab_mb.athinput/he_slab_mb_$k.athinput}
  ./run.sh $n new $x "$@"; }
gt tb3_Fx2r4 Fx2r4 tbd $TB3 & gt tb3_Ux1r4 Ux1r4 tbd $TB3 & gt tb3_Fx12r4 Fx12r4 tbd $TB3 & wait
gt tb3_T3r4 T3r4 tb $TB3 & gt tb3_Uoddr3 Uoddr3 tbd $TB3 & gt tb3_T3r2 T3r2 tb $TB3 & gt tb3_T3r1 T3r1 tb $TB3 & gt tb3_Ux2r1 Ux2r1 tbd $TB3 & wait
gt tb3_Ux2b2r2 Ux2b2r2 tbd $TB3 & gt tb2_T3r4 T3r4 tb rad_m1/vet_mb_tblock=2 & gt tb2_Ux2b2r2 Ux2b2r2 tb rad_m1/vet_mb_tblock=2 & gt tb3_U1 U1 tbd $TB3 & wait
gt tb3t_T3r4 T3r4 tb $TB3 rad_m1/vet_mb_tblock_tj=3 rad_m1/vet_mb_tblock_tk=5 & gt tb3t_T3r1 T3r1 tb $TB3 rad_m1/vet_mb_tblock_tj=4 rad_m1/vet_mb_tblock_tk=3 & gt tb3t_Fx2r4 Fx2r4 tbd $TB3 rad_m1/vet_mb_tblock_tj=2 & gt tb3no_T3r4 T3r4 tb $TB3 rad_m1/vet_mb_tblock_ovl=false & wait
gt tb3na_T3r4 T3r4 tb $TB3 rad_m1/vet_mb_agg=false & gt tb4_T3r4 T3r4 tb rad_m1/vet_mb_tblock=4 rad_m1/vet_mb_halo=4 & gt tb4d_Fx12r4 Fx12r4 tbd rad_m1/vet_mb_tblock=4 & wait
gt tb4d_T3r4 T3r4 tbd rad_m1/vet_mb_tblock=4 & gt tb3h1_T3r4 T3r4 tb $TB3 rad_m1/vet_mb_halo=1 & gt tb3t1_Fx12r4 Fx12r4 tbd $TB3 rad_m1/vet_mb_tblock_team=1 rad_m1/vet_mb_tblock_tj=3 & wait
./run.sh hy2r_T3r4 ref ${A[T3r4]/he3d.athinput/he3d_tb.athinput} rad_m1/vet_mb_agroup=2 & gt hy2t_T3r4 T3r4 tb rad_m1/vet_mb_agroup=2 $TB3 rad_m1/vet_mb_tblock_tj=5 & ./run.sh hy2r_Fx12r4 ref ${A[Fx12r4]/he_slab_mb.athinput/he_slab_mb_tbd.athinput} rad_m1/vet_mb_agroup=2 & wait
# bases: vet_mb_tblock = 1 on the same inputs (the bin/rst headers then differ only in
# the echoed override values)
gt b1_T3r4 T3r4 tb & gt b1_T3r2 T3r2 tb & gt b1_T3r1 T3r1 tb & gt b1_Ux2b2r2 Ux2b2r2 tb & wait
gt b1d_Fx2r4 Fx2r4 tbd & gt b1d_Ux1r4 Ux1r4 tbd & gt b1d_Fx12r4 Fx12r4 tbd & wait
gt b1d_Uoddr3 Uoddr3 tbd & gt b1d_Ux2r1 Ux2r1 tbd & gt b1d_Ux2b2r2 Ux2b2r2 tbd & gt b1d_U1 U1 tbd & gt b1d_T3r4 T3r4 tbd & wait
gt hy2t_Fx12r4 Fx12r4 tbd rad_m1/vet_mb_agroup=2 $TB3 & wait
echo "== A. default arms (tblock off), new vs ref: every file"
for a in U1 Fx2r4 Ux2r1 Ux2b2r2 Ux1r4 Fx12r4 Uoddr3 T3r4 T3r2 T3r1; do ./compare.sh ref_$a new_$a; done
echo "== B. tblock arms: vs ref (hst), vs the tblock = 1 base on the same input (all bytes)"
for x in tb3_Fx2r4:Fx2r4:d tb3_Ux1r4:Ux1r4:d tb3_Fx12r4:Fx12r4:d tb3_T3r4:T3r4: tb3_Uoddr3:Uoddr3:d tb3_T3r2:T3r2: tb3_T3r1:T3r1: tb3_Ux2r1:Ux2r1:d tb3_Ux2b2r2:Ux2b2r2:d tb3_U1:U1:d tb2_T3r4:T3r4: tb2_Ux2b2r2:Ux2b2r2: tb3t_T3r4:T3r4: tb3t_T3r1:T3r1: tb3t_Fx2r4:Fx2r4:d tb3no_T3r4:T3r4: tb3na_T3r4:T3r4: tb4_T3r4:T3r4: tb4d_Fx12r4:Fx12r4:d tb4d_T3r4:T3r4:d tb3h1_T3r4:T3r4: tb3t1_Fx12r4:Fx12r4:d; do
  IFS=: read n a k <<< "$x"; ./compare.sh ref_$a $n; ./compare.sh b1${k}_$a $n; done
echo "== B0. bases (tblock = 1 spelled out) vs ref: hst"
for x in b1_T3r4:T3r4 b1_T3r2:T3r2 b1_T3r1:T3r1 b1_Ux2b2r2:Ux2b2r2 b1d_Fx2r4:Fx2r4 b1d_Ux1r4:Ux1r4 b1d_Fx12r4:Fx12r4 b1d_Uoddr3:Uoddr3 b1d_Ux2r1:Ux2r1 b1d_Ux2b2r2:Ux2b2r2 b1d_U1:U1 b1d_T3r4:T3r4; do ./compare.sh ref_${x#*:} ${x%%:*}; done
echo "== C. hybrid G=2 with tblock vs hybrid G=2 HEAD"
./compare.sh hy2r_T3r4 hy2t_T3r4; ./compare.sh hy2r_Fx12r4 hy2t_Fx12r4
