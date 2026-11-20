#!/bin/bash
# apply_v202.sh (INTEG6): the USER's v20.2 key decision (2026-11-20 ~13:20 PST), run ONLY at the v20.2 push. Backups *.bak_v201.
#   hydro/floor_pay=rad_column for ALL stars (HG, RSG, R136): the hydro C2P floor energy is paid from the radiation E
#   (ELEAK; closes the rule-1 energy residual; result-changing). Added to uniform/keys/common_scheme.txt; otherwise = v20.
#   RAVEN=1 also writes the Raven speedhead keys (from integ/keys_v202 on Raven).
S=/orion/ptmp/jinma/clean_1115/stag2; U=$S/uniform/keys
for s in hg rsg r136; do
  f=$S/${s}_latest.txt; cp $f $f.bak_v201
  cat $U/common_scheme.txt $U/${s}_problem.txt | grep -v '^#' | grep -v '^[[:space:]]*$' > $f
  a=$S/integ/keys/k1_${s}.args; cp $a $a.bak_v201; tr '\n' ' ' < $f > $a
done
if [ "${RAVEN:-0}" = 1 ]; then
  ssh -S ~/.ssh/raven.sock raven01i "cd /raven/ptmp/jinma/stag2/speedhead; for s in hg rsg r136; do cp keys_\$s.txt keys_\$s.txt.bak_v201; cp /raven/ptmp/jinma/stag2/integ/keys_v202/\${s}_v202.txt keys_\$s.txt; done"
fi
for s in hg rsg r136; do diff $S/${s}_latest.txt.bak_v201 $S/${s}_latest.txt; done; echo "apply_v202 done"
