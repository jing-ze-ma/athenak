#!/bin/bash
# apply_v205.sh (INTEG9): keys of v20.5 for HG / RSG / R136 = v20.4 option B (USER 2026-11-20 ~20:45 PST: the v20.3 behaviour
# named; result-changing keys enter production only after their own long runs) + rad_m1/clean_stag_flim=false named (flim is a
# named key, default off, USER 2026-11-21: the production restarts keep it off; naming it follows the UNIFORM rule "name every
# scheme value"; the fresh 1121b production decks name clean_stag_flim=true themselves). clean_stag_thinf_held (BEAM3 C1) is not named (default off; needs THINF, off here).
# Run ONLY at the push (it supersedes apply_v204*.sh; backups *.bak_v203). DRY=<dir>: write the three keys files to <dir> only.
# RAVEN=1 also writes the Raven speedhead keys: its v20.3 keys_<star>.txt + the same NEW keys (the v20.4 Raven keys_v204/ are
# option A and must NOT be copied).
S=/orion/ptmp/jinma/clean_1115/stag2; U=$S/uniform/keys
NEW="rad_m1/clean_stag_thinf=false rad_m1/clean_maxit=200 rad_m1/clean_conv_max=1e-4 rad_m1/clean_maxnl=200 rad_m1/clean_lin_pass=true rad_m1/clean_pc_sweep=rb rad_m1/clean_stag_flim=false"
gen() {  # gen OUTDIR: common + NEW + <star>_problem -> OUTDIR/<star>_latest.txt (first value of a key wins)
  for s in hg rsg r136; do
    { grep -v '^#' $U/common_scheme.txt; printf '%s\n' $NEW; grep -v '^#' $U/${s}_problem.txt; } \
      | grep -v '^[[:space:]]*$' | awk -F= '!seen[$1]++' > $1/${s}_latest.txt; done; }
if [ -n "$DRY" ]; then mkdir -p $DRY; gen $DRY; for s in hg rsg r136; do echo "== $s"; diff <(tr ' ' '\n' < $S/${s}_latest.txt | sort) <(sort $DRY/${s}_latest.txt); done; exit 0; fi
cp $U/common_scheme.txt $U/common_scheme.txt.bak_v203
{ echo "# --- v20.5 (= v20.4 option B + flim named off): the v20.3 behaviour named for the production restarts (clean_maxit 800 / clean_conv_max 1e-5 are plain defaults and would change restarts; THINF, pc_sweep sym, clean_stag_flim stay off) ---"; printf '%s\n' $NEW; } >> $U/common_scheme.txt
for s in hg rsg r136; do
  f=$S/${s}_latest.txt; cp $f $f.bak_v203
  cat $U/common_scheme.txt $U/${s}_problem.txt | grep -v '^#' | grep -v '^[[:space:]]*$' | awk -F= '!seen[$1]++' > $f
  a=$S/integ/keys/k1_${s}.args; cp $a $a.bak_v203; tr '\n' ' ' < $f > $a
done
if [ "${RAVEN:-0}" = 1 ]; then
  ssh -S ~/.ssh/raven.sock raven01i "cd /raven/ptmp/jinma/stag2/speedhead; for s in hg rsg r136; do cp keys_\$s.txt keys_\$s.txt.bak_v203; { cat keys_\$s.txt.bak_v203; printf '%s\n' $NEW; } | awk -F= '!seen[\$1]++' > keys_\$s.txt; done"
fi
for s in hg rsg r136; do diff $S/${s}_latest.txt.bak_v203 $S/${s}_latest.txt; done; echo "apply_v205 done"
