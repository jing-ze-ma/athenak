#!/bin/bash
# apply_v203.sh (INTEG7): the USER's v20.3 key decision (2026-11-20 ~14:55 PST), run ONLY at the v20.3 push. Backups *.bak_v202.
#   rad_m1/clean_cs_lz_tel=true for HG and RSG (cubed sphere, staggered; LZTEL3 fix 74c8a8bb1: the net lateral-stress
#   torque defect goes to the gas as a least-KE solid-body spin; result-changing). Added to uniform/keys/{hg,rsg}_problem.txt.
#   rad_m1/vet_rotate_avg_ema=128 + rad_m1/vet_rotate_avg_ema_change_form=step for HG and RSG (USER 2026-11-20): what the
#   production restarts hegiant.00184 / rsg.01510 carry in their files; named so a restart from any file runs one form.
#   R136 unchanged (its restart r136.00097 carries no EMA: K4 window). vet_rotate_avg_ema_wb stays unnamed (off; USER: a Caltech wb on/off pair decides later).
#   RAVEN=1 also writes the Raven speedhead keys (from integ/keys_v203 on Raven).
S=/orion/ptmp/jinma/clean_1115/stag2; U=$S/uniform/keys
for s in hg rsg; do p=$U/${s}_problem.txt; cp $p $p.bak_v202
  for kv in rad_m1/vet_rotate_avg_ema=128 rad_m1/vet_rotate_avg_ema_change_form=step rad_m1/clean_cs_lz_tel=true; do
    grep -q "^${kv%%=*}=" $p || echo $kv >> $p; done; done
for s in hg rsg r136; do
  f=$S/${s}_latest.txt; cp $f $f.bak_v202
  cat $U/common_scheme.txt $U/${s}_problem.txt | grep -v '^#' | grep -v '^[[:space:]]*$' > $f
  a=$S/integ/keys/k1_${s}.args; cp $a $a.bak_v202; tr '\n' ' ' < $f > $a
done
if [ "${RAVEN:-0}" = 1 ]; then
  ssh -S ~/.ssh/raven.sock raven01i "cd /raven/ptmp/jinma/stag2/speedhead; for s in hg rsg r136; do cp keys_\$s.txt keys_\$s.txt.bak_v202; cp /raven/ptmp/jinma/stag2/integ/keys_v203/\${s}_v203.txt keys_\$s.txt; done"
fi
for s in hg rsg r136; do diff $S/${s}_latest.txt.bak_v202 $S/${s}_latest.txt; done; echo "apply_v203 done"
