#!/bin/bash
# switch_k9.sh RUN NODE HOLDER [CFL]   PREPARED, NOT RUN (orion ki / kj): HOLD until the USER's yes is relayed by the coordinator.
# Production switch of RUN to B_opac3 (athena_opac3_772ee5ac6, md5 3541666326ecc8b62fa14547f5ffd5b2) + clean_conv_max 1e-4 + clean_face_opac harmonic + clean_stag_opac_update + cfl (0.9; fallback 0.6).
# B_j3 ignores STOP / SIGUSR1, so the switch point is the next periodic rst: wait until a rst NEWER than the current one is complete (size stable over 60 s),
# then cancel the production step (scancel of the p_<RUN> step only, never a holder), check nvidia-smi is empty, and restart from that rst via prodg_k9.sh (log_k9.<n>, chunks.txt with md5s).
J=/resnick/groups/carnegie_poc/jingze; P=$J/prod_fresh_1121b; R=$1; NODE=$2; HOLDER=$3; CFL=${4:-0.9}; B=$J/wt_orion1121/bin/athena_opac3_772ee5ac6; D=$P/$R
[ "$(md5sum < $B | cut -c1-32)" = 3541666326ecc8b62fa14547f5ffd5b2 ] || { echo "B_opac3 md5 mismatch"; exit 1; }
cur=$(ls -t $D/rst/*.rst | head -1); echo "$(date '+%F %T') waiting for the next rst after $(basename $cur)"
until new=$(ls -t $D/rst/*.rst | head -1); [ "$new" != "$cur" ] && s1=$(stat -c %s $new) && sleep 60 && [ "$s1" = "$(stat -c %s $new)" ]; do sleep 30; done
st=$(squeue -s -h -o '%i %j' -u $USER | awk -v j=p_$R '$2==j{print $1}'); echo "$(date '+%F %T') rst $(basename $new) complete ($(stat -c %s $new) B); cancelling step $st"; [ -n "$st" ] && scancel $st
until [ $(srun --overlap --jobid=$HOLDER -w $NODE -N1 -n1 -c1 -t 2 -J p_smi nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -cE '^[0-9]+') = 0 ]; do sleep 10; done
while IFS= read -r l; do export "$l"; done < $P/i1_m35mms_f1k.env 2>/dev/null
export NODE HOLDER KEEPRST=3 END="2026-11-21 17:25" RST_K9=$new
exec $P/prodg_k9.sh $R $B $(awk '/DECK=/{for(i=1;i<=NF;i++) if($i ~ /^DECK=/){sub("DECK=","",$i); print $i}}' $D/chunks.txt | head -1) rad_m1/clean_conv_max=1.0e-4 rad_m1/clean_face_opac=harmonic rad_m1/clean_stag_opac_update=true time/cfl_number=$CFL
