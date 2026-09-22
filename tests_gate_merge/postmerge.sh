#!/bin/bash -l
# post-merge gate: rebuild box_convection at HEAD, box G1 (modes 3, 0) vs the references,
# then the 1-D He column with implicit_x1 for 300 s as a smoke test of the merged tree
module purge; module load gcc/14 cmake/4.0
A=/viper/u2/jinma/ATHENAK/athenak; G=$A/tests_gate_merge
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
make -C $A/build_cpu_box -j 12 > $G/postmerge_build.log 2>&1 || { echo BUILD_FAILED; exit 1; }
for M in 3 0; do
  D=$G/g1_post_m$M; rm -rf $D; mkdir -p $D; cd $D
  OPT=""; [ "$M" = "0" ] && OPT="problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false"
  $A/build_cpu_box/src/athena -i $IN mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 \
    meshblock/nx2=8 meshblock/nx3=8 time/nlim=50 $OPT > run.log 2>&1
  for f in feczrt.hydro.hst feczrt.user.hst column_used.txt rt_surface.bin rt_profile.bin; do
    cmp -s $f $G/g1_new_m$M/$f && echo "m$M $f IDENTICAL" || echo "m$M $f DIFFER"
  done
  rm -rf bin cbin* rst
done
D=$G/hecol_impl; rm -rf $D; mkdir -p $D; cd $D
cp $A/tests_m1/runs_2a/m1_rad_ic_V2.txt . 2>/dev/null
$A/build_cpu_box/src/athena -i $A/tests_m1/runs_3a/he_box_m1_1d_impl.athinput time/tlim=300 \
  > run.log 2>&1; echo "hecol exit=$?"
tail -3 run.log; ls *.user.hst >/dev/null 2>&1 && tail -1 *.user.hst | cut -c1-200
rm -rf bin cbin* rst
echo POSTMERGE_DONE
