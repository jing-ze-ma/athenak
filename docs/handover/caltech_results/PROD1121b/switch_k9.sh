#!/bin/bash
# switch_k9.sh RUN NODE HOLDER   PREPARED, NOT RUN (orion ki / kj / kl): HOLD until the USER's yes is relayed by the coordinator.
# Production switch of RUN to B_fdef (caltech-test-fdef 03b027502, athena_fdef_03b027502; md5 checked against its .md5 at launch) with
#   rst keys + clean_conv_max 1e-4 + clean_face_opac harmonic + clean_stag_opac_update true + clean_stag_flim false + clean_stag_fwb true + clean_stag_fcone true + cfl 0.9 (no damp, no gate).
# The current binary (B_j3) ignores STOP / SIGUSR1, so the switch point is the next periodic rst: wait until a rst NEWER than the current one is complete (size stable 60 s),
# cancel only the p_<RUN> step (never a holder), check nvidia-smi is empty, restart from that rst via prodg_k9.sh (log_k9.<n>, chunks.txt with md5s, KEEPRST 3, END 17:25).
# Fallback ladder (orion kl), from the SAME rst, on any CW / FATAL / DT_RETRY in the first 100 cycles: (1) cfl 0.6, same keys; (2) cfl 0.6 with flim true, fwb false, fcone false.
J=/resnick/groups/carnegie_poc/jingze; P=$J/prod_fresh_1121b; R=$1; NODE=$2; HOLDER=$3; B=$J/wt_orion1121/bin/athena_fdef_03b027502; D=$P/$R
[ -s $B.md5 ] && [ "$(md5sum < $B | cut -c1-32)" = "$(cut -c1-32 $B.md5)" ] || { echo "B_fdef missing or md5 mismatch"; exit 1; }
DECK=$(awk '/DECK=/{for(i=1;i<=NF;i++) if($i ~ /^DECK=/){sub("DECK=","",$i); print $i}}' $D/chunks.txt | head -1)
K0="rad_m1/clean_conv_max=1.0e-4 rad_m1/clean_face_opac=harmonic rad_m1/clean_stag_opac_update=true"
L1="$K0 rad_m1/clean_stag_flim=false rad_m1/clean_stag_fwb=true rad_m1/clean_stag_fcone=true time/cfl_number=0.9"
L2="$K0 rad_m1/clean_stag_flim=false rad_m1/clean_stag_fwb=true rad_m1/clean_stag_fcone=true time/cfl_number=0.6"
L3="$K0 rad_m1/clean_stag_flim=true rad_m1/clean_stag_fwb=false rad_m1/clean_stag_fcone=false time/cfl_number=0.6"
cur=$(ls -t $D/rst/*.rst | head -1); echo "$(date '+%F %T') waiting for the next rst after $(basename $cur)"
until new=$(ls -t $D/rst/*.rst | head -1); [ "$new" != "$cur" ] && s1=$(stat -c %s $new) && sleep 60 && [ "$s1" = "$(stat -c %s $new)" ]; do sleep 30; done
cp -p $new $D/switch_k9_from.$(basename $new) 2>/dev/null || true   # keep the switch rst safe from KEEPRST pruning
stopstep() { local st=$(squeue -s -h -o '%i %j' -u $USER | awk -v j=p_$R '$2==j{print $1}'); [ -n "$st" ] && scancel $st; until [ $(srun --overlap --jobid=$HOLDER -w $NODE -N1 -n1 -c1 -t 2 -J p_smi nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -cE '^[0-9]+') = 0 ]; do sleep 10; done; }
echo "$(date '+%F %T') rst $(basename $new) complete ($(stat -c %s $new) B); stopping the production step"; stopstep
while IFS= read -r l; do export "$l"; done < $P/i1_m35mms_f1k.env 2>/dev/null
export NODE HOLDER KEEPRST=3 END="2026-11-21 17:25" RST_K9=$D/switch_k9_from.$(basename $new)
for lvl in 1 2 3; do eval KEYS=\$L$lvl; echo "$(date '+%F %T') launch level $lvl: $KEYS"
  setsid $P/prodg_k9.sh $R $B $DECK $KEYS < /dev/null > $D/switch_k9_l$lvl.out 2>&1 & PG=$!
  sleep 60; LG=$(ls -t $D/log_k9.* 2>/dev/null | head -1); c0=""; bad=""
  while kill -0 $PG 2>/dev/null; do c=$(grep -oE 'cycle=[0-9]+' $LG | tail -1 | cut -d= -f2); [ -z "$c0" ] && c0=$(grep -m1 -oE 'cycle=[0-9]+' $LG | cut -d= -f2)
    grep -qE 'CLEAN WARNING|FATAL|DT_RETRY' $LG && { bad=$(grep -m1 -E 'CLEAN WARNING|FATAL|DT_RETRY' $LG | cut -c1-160); break; }
    [ -n "$c0" ] && [ -n "$c" ] && [ $((c - c0)) -ge 100 ] && break; sleep 20; done
  if [ -z "$bad" ]; then echo "$(date '+%F %T') level $lvl clean over the first 100 cycles; production continues on it"; exit 0; fi
  echo "$(date '+%F %T') level $lvl: '$bad' in the first 100 cycles -> next level"; stopstep; done
echo "$(date '+%F %T') all three levels hit CW / FATAL / DT_RETRY in the first 100 cycles; production of $R is STOPPED, report to the coordinator"
