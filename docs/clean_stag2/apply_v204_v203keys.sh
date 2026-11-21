#!/bin/bash
# apply_v204_v203keys.sh (INTEG8): the v20.4 keys with the v20.3 BEHAVIOUR named (USER option B; CM5_hg FATAL with option A), run ONLY at the v20.4 push (backups *.bak_v203). DRY=<dir>: write the three keys files
# to <dir> only (nothing else touched).
#   v20.4 code defaults (USER 2026-11-20 ~15:55 PST), named in the keys anyway (UNIFORM: name every scheme value):
#     rad_m1/clean_stag_thinf=true + rad_m1/clean_maxit=800 (beam fix THINF v7 G=1; result change for HG/RSG/R136: the
#     production keys did not name it), rad_m1/vet_rotate_avg_ema=128 + _change_form=step
#     (EMA moved from hg/rsg_problem.txt to common_scheme.txt: R136 takes it too, USER 2026-11-20: result change for R136);
#     rad_m1/clean_cs_lz_tel=true stays in hg/rsg_problem.txt (cubed sphere only; R136 is spherical polar: FATAL there; and
#     its code default is off, USER 2026-11-20 ~17:25 PST).
#   rad_m1/clean_maxnl=200 + rad_m1/clean_lin_pass=true (USER 2026-11-20: the old production values, no longer inherited from
#     implicit_maxit / implicit_lin_maxit; bitwise for the production restarts).
#   rad_m1/clean_conv_max=1e-5 (USER 2026-11-20 ~16:22 PST: production convergence ceiling; the restarts' implicit_tol_rel_max=1e-4
#     fed it before; result change for all three).
#   RAVEN=1 also writes the Raven speedhead keys (from integ/keys_v204 on Raven).
#   OPTION B: every production star keeps its v20.3 results on the v20.4 binary: clean_stag_thinf=false, clean_maxit=200,
#   clean_conv_max=1e-4 (both were the restarts' implicit_lin_maxit / implicit_tol_rel_max), clean_pc_sweep=rb, clean_maxnl=200,
#   clean_lin_pass=true named; EMA 128/step and clean_cs_lz_tel stay in hg/rsg_problem.txt, R136 none (as v20.3). Verified: V204O
#   (v20.3 keys + maxit 200 + conv 1e-4) == V203 bitwise on RSG (Caltech) and R136 (Raven).
S=/orion/ptmp/jinma/clean_1115/stag2; U=$S/uniform/keys; NEW="rad_m1/clean_stag_thinf=false rad_m1/clean_maxit=200 rad_m1/clean_conv_max=1e-4 rad_m1/clean_maxnl=200 rad_m1/clean_lin_pass=true rad_m1/clean_pc_sweep=rb"
MV="rad_m1/__none__="
gen() {  # gen OUTDIR: common (+ the v20.4 block) and <star>_problem (without the moved keys) -> OUTDIR/<star>_latest.txt
  for s in hg rsg r136; do
    { grep -v '^#' $U/common_scheme.txt; printf '%s\n' $NEW; grep -v '^#' $U/${s}_problem.txt | grep -Ev "^($MV)"; } \
      | grep -v '^[[:space:]]*$' | awk -F= '!seen[$1]++' > $1/${s}_latest.txt; done; }
if [ -n "$DRY" ]; then mkdir -p $DRY; gen $DRY; for s in hg rsg r136; do echo "== $s"; diff <(tr ' ' '\n' < $S/${s}_latest.txt | sort) <(sort $DRY/${s}_latest.txt); done; exit 0; fi
cp $U/common_scheme.txt $U/common_scheme.txt.bak_v203
{ echo "# --- v20.4 option B: the v20.3 behaviour named (the v20.4 code defaults clean_maxit 800 / clean_conv_max 1e-5 are plain defaults and would change restarts; THINF and pc_sweep sym stay off) ---"; printf '%s\n' $NEW; } >> $U/common_scheme.txt
for s in hg rsg; do cp $U/${s}_problem.txt $U/${s}_problem.txt.bak_v203; grep -Ev "^($MV)" $U/${s}_problem.txt.bak_v203 > $U/${s}_problem.txt; done
for s in hg rsg r136; do
  f=$S/${s}_latest.txt; cp $f $f.bak_v203
  cat $U/common_scheme.txt $U/${s}_problem.txt | grep -v '^#' | grep -v '^[[:space:]]*$' > $f
  a=$S/integ/keys/k1_${s}.args; cp $a $a.bak_v203; tr '\n' ' ' < $f > $a
done
if [ "${RAVEN:-0}" = 1 ]; then
  ssh -S ~/.ssh/raven.sock raven01i "cd /raven/ptmp/jinma/stag2/speedhead; for s in hg rsg r136; do cp keys_\$s.txt keys_\$s.txt.bak_v203; cp /raven/ptmp/jinma/stag2/integ/keys_v204/\${s}_v204.txt keys_\$s.txt; done"
fi
for s in hg rsg r136; do diff $S/${s}_latest.txt.bak_v203 $S/${s}_latest.txt; done; echo "apply_v204_v203keys done"
