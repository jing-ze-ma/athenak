#!/bin/bash -l
# gates.sh <tag>  -- run all bitwise gates with build_clean_* binaries into tests_cleanup_0922/<tag>
module purge >/dev/null 2>&1; module load gcc/14 cmake/4.0 >/dev/null 2>&1
A=/viper/u2/jinma/ATHENAK/athenak; G=$A/tests_cleanup_0922; T=$G/$1
rm -rf $T; mkdir -p $T
BOX=$A/build_clean_box_convection/src/athena
RG=$A/build_clean_red_giant/src/athena
DHJ=$A/build_clean_deep_hot_jupiter_rt/src/athena
CS=$A/build_clean_cs_test/src/athena
SP=$A/build_clean_sp_test/src/athena

# (a1) M1 He slab, box_convection pgen, tlim 3
D=$T/a1_m1slab; mkdir -p $D; cd $D
cp $A/tests_m1/runs_3d_edd/m1_rad_ic_V3edd.txt .
timeout 900 $BOX -i $A/tests_m1/runs_3d_edd/he_slab_m1_2d_V3edd.athinput time/tlim=3.0 \
  > run.log 2>&1; echo "a1 exit=$?"

# (a2) box G1 configuration from tests_gate_merge/postmerge.sh (modes 3 and 0)
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
for M in 3 0; do
  D=$T/a2_g1_m$M; mkdir -p $D; cd $D
  OPT=""; [ "$M" = "0" ] && OPT="problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false"
  timeout 900 $BOX -i $IN mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 \
    meshblock/nx2=8 meshblock/nx3=8 time/nlim=50 $OPT > run.log 2>&1; echo "a2_m$M exit=$?"
done

# (b1) red giant 1-D radiative column
D=$T/b1_rg_col; mkdir -p $D; cd $D
timeout 900 $RG -i $A/inputs/hydro/red_giant_column.athinput time/nlim=300 \
  problem/opac_table=$A/data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt > run.log 2>&1
echo "b1 exit=$?"
# (b2) red giant 3-D production input, reduced angular resolution
D=$T/b2_rg_fofc; mkdir -p $D; cd $D
timeout 900 $RG -i $A/inputs/hydro/red_giant_fofc.athinput time/nlim=20 \
  mesh/nx2=8 mesh/nx3=8 meshblock/nx2=8 meshblock/nx3=8 > run.log 2>&1
echo "b2 exit=$?"

# (c) dhj ck spherical, 20 cycles, both ck_spherical values
for S in false true; do
  D=$T/c_dhj_$S; mkdir -p $D; cd $D
  timeout 900 $DHJ -i $A/inputs/tests/dhj_ck_spherical.athinput time/nlim=20 \
    problem/ck_spherical=$S > run.log 2>&1; echo "c_$S exit=$?"
done

# (d) cubed sphere blast (hydro + mhd)
D=$T/d_cs_blast; mkdir -p $D; cd $D
timeout 900 $CS -i $A/inputs/tests/cubed_sphere_blast.athinput time/nlim=20 > run.log 2>&1
echo "d1 exit=$?"
D=$T/d_cs_mhd_blast; mkdir -p $D; cd $D
timeout 900 $CS -i $A/inputs/tests/cubed_sphere_mhd_blast.athinput time/nlim=20 > run.log 2>&1
echo "d2 exit=$?"

# (e) spherical polar MHD blast (polar treatment)
D=$T/e_sp_blast_mhd; mkdir -p $D; cd $D
timeout 900 $SP -i $A/inputs/tests/spherical_polar_blast_mhd.athinput time/nlim=20 > run.log 2>&1
echo "e exit=$?"

# checksums of every payload except logs
# NOTE: rst/ is excluded -- restart files embed the effective parameter list, which
# legitimately changes when a parameter is deleted; the data payload is covered by bin/.
cd $T && find . -type f ! -name 'run.log' ! -name '*.txt' ! -path './*/rst/*' -print0 | sort -z \
  | xargs -0 md5sum > $G/$1.md5
echo "GATES_DONE $1  ($(wc -l < $G/$1.md5) files)"
