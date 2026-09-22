#!/bin/bash -l
# gates_old.sh <tag> -- gates.sh, but with every default flipped on 2026-09-22 written
# back into the input file EXPLICITLY at its OLD value.  Used for gate (a): the new
# binary must then be bitwise the reference.  Only inputs that were SILENT about a
# flipped switch are patched; where the input already states the switch (he_box_w8's
# rt_col3_ex_iter = true, dhj_ck_spherical's ck_spherical / ck_beam_sph / rt_split) the
# old effective value is what the file already says.
module purge >/dev/null 2>&1; module load gcc/14 cmake/4.0 >/dev/null 2>&1
A=/viper/u2/jinma/ATHENAK/athenak; G=$A/tests_cleanup_0922; T=$G/$1
rm -rf $T; mkdir -p $T/inputs
BOX=$A/build_clean_box_convection/src/athena
RG=$A/build_clean_red_giant/src/athena
DHJ=$A/build_clean_deep_hot_jupiter_rt/src/athena
CS=$A/build_clean_cs_test/src/athena
SP=$A/build_clean_sp_test/src/athena

# patch <src> <dst> <line>... : insert the lines just after the "<problem>" header
patch_in() {
  python3 - "$@" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().split('\n')
i = [k for k, l in enumerate(lines) if l.startswith('<problem>')][0]
for n, extra in enumerate(sys.argv[3:]):
    lines.insert(i + 1 + n, extra)
open(dst, 'w').write('\n'.join(lines))
PY
}
M1=$A/tests_m1/runs_3d_edd/he_slab_m1_2d_V3edd.athinput
patch_in $M1 $T/inputs/m1slab.athinput "rt_col3_ex_iter = false"
patch_in $A/inputs/hydro/red_giant_column.athinput $T/inputs/rg_col.athinput \
  "rt_col3_ex_iter = false" "mlt_split_deposit = false"
patch_in $A/inputs/hydro/red_giant_fofc.athinput $T/inputs/rg_fofc.athinput \
  "rt_col3_ex_iter = false" "mlt_split_deposit = false"
patch_in $A/inputs/tests/dhj_ck_spherical.athinput $T/inputs/dhj.athinput \
  "rt_use_cons = false"

D=$T/a1_m1slab; mkdir -p $D; cd $D
cp $A/tests_m1/runs_3d_edd/m1_rad_ic_V3edd.txt .
timeout 900 $BOX -i $T/inputs/m1slab.athinput time/tlim=3.0 > run.log 2>&1
echo "a1 exit=$?"

IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
for M in 3 0; do
  D=$T/a2_g1_m$M; mkdir -p $D; cd $D
  OPT=""; [ "$M" = "0" ] && OPT="problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false"
  timeout 900 $BOX -i $IN mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 \
    meshblock/nx2=8 meshblock/nx3=8 time/nlim=50 $OPT > run.log 2>&1; echo "a2_m$M exit=$?"
done

D=$T/b1_rg_col; mkdir -p $D; cd $D
timeout 900 $RG -i $T/inputs/rg_col.athinput time/nlim=300 \
  problem/opac_table=$A/data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt > run.log 2>&1
echo "b1 exit=$?"
D=$T/b2_rg_fofc; mkdir -p $D; cd $D
timeout 900 $RG -i $T/inputs/rg_fofc.athinput time/nlim=20 \
  mesh/nx2=8 mesh/nx3=8 meshblock/nx2=8 meshblock/nx3=8 > run.log 2>&1
echo "b2 exit=$?"

for S in false true; do
  D=$T/c_dhj_$S; mkdir -p $D; cd $D
  timeout 900 $DHJ -i $T/inputs/dhj.athinput time/nlim=20 \
    problem/ck_spherical=$S > run.log 2>&1; echo "c_$S exit=$?"
done

D=$T/d_cs_blast; mkdir -p $D; cd $D
timeout 900 $CS -i $A/inputs/tests/cubed_sphere_blast.athinput time/nlim=20 > run.log 2>&1
echo "d1 exit=$?"
D=$T/d_cs_mhd_blast; mkdir -p $D; cd $D
timeout 900 $CS -i $A/inputs/tests/cubed_sphere_mhd_blast.athinput time/nlim=20 \
  > run.log 2>&1
echo "d2 exit=$?"

D=$T/e_sp_blast_mhd; mkdir -p $D; cd $D
timeout 900 $SP -i $A/inputs/tests/spherical_polar_blast_mhd.athinput time/nlim=20 \
  > run.log 2>&1
echo "e exit=$?"

cd $T && find . -type f ! -name 'run.log' ! -name '*.txt' ! -path './*/rst/*' \
  ! -path './inputs/*' -print0 | sort -z | xargs -0 md5sum > $G/$1.md5
echo "GATES_DONE $1  ($(wc -l < $G/$1.md5) files)"
