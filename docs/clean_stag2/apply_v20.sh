#!/bin/bash
# apply_v20.sh (INTEG3): the USER's v20 key decision (2026-11-20), run ONLY at the v20 cut. Backups *.bak_v19.
#   UNIFORM common scheme set (stag2/uniform/UNIFORM.md): every LATEST keys file = uniform/keys/common_scheme.txt +
#   uniform/keys/<star>_problem.txt (comments dropped). Every scheme value is named explicitly (restarts keep the old
#   defaults otherwise), incl. clean_conv_kappa=1.0 (= the v20 code default), time2_k1_stiff=0, time/restart_refill_ghosts=false,
#   clean_rot_fuse=true for all; the 5 inert keys (time2_adm_poslim, time2_dtr_poslim, zalesak_cap_rounds, vet_col_lat_posfix,
#   time2_vet_x2_fields; preflight "no effect" since 49200b4c1) and the other inert keys are retired.
#   Hardware keys (uniform/keys/hardware_rules.txt, e.g. time/dt_retry_host=m1 for HG 12 ranks / 3 A100 nodes) are NOT in the
#   scheme files: name them in the job's hardware deck.
#   time2_k1_stiff=0 for ALL stars (USER 2026-11-20: HG ~7 % cost over 30 cycles accepted; RSG +84 % with 1.0).
#   RAVEN=1 also writes the Raven speedhead keys (Raven table paths, integ/keys_v20 on Raven).
S=/orion/ptmp/jinma/clean_1115/stag2; U=$S/uniform/keys
for s in hg rsg r136; do
  f=$S/${s}_latest.txt; cp $f $f.bak_v19
  cat $U/common_scheme.txt $U/${s}_problem.txt | grep -v '^#' | grep -v '^[[:space:]]*$' > $f
  a=$S/integ/keys/k1_${s}.args; cp $a $a.bak_v19; tr '\n' ' ' < $f > $a
done
if [ "${RAVEN:-0}" = 1 ]; then
  ssh -S ~/.ssh/raven.sock raven01i "cd /raven/ptmp/jinma/stag2/speedhead; for s in hg rsg r136; do cp keys_\$s.txt keys_\$s.txt.bak_v19; cp /raven/ptmp/jinma/stag2/integ/keys_v20/\${s}_v20.txt keys_\$s.txt; done"
fi
for k in time2_adm_poslim time2_dtr_poslim zalesak_cap_rounds vet_col_lat_posfix time2_vet_x2_fields; do grep -H "/$k=" $S/{hg,rsg,r136}_latest.txt; done
grep -H "clean_conv_kappa\|time2_k1_stiff\|restart_refill_ghosts\|clean_rot_fuse" $S/{hg,rsg,r136}_latest.txt; echo "apply_v20 done"
