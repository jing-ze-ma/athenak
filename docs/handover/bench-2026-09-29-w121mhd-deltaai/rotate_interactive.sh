#!/bin/bash
# rotate_interactive.sh: move the scan arms one at a time from ghx4 to ghx4-interactive (qos_ghx4int allows 1 job
# queued or running per user; user 09-29: one at a time, not packed). Runs on the login node.
# rotation.txt: "ARM GHX4_JOBID NRANK LIMIT", shortest first. When no interactive job of ours is queued/running,
# take the next arm whose ghx4 job is still PENDING, cancel that job and submit the arm to ghx4-interactive.
# An arm whose ghx4 job already started is left there. Exits when no pending arm is left. Log: rotation.log
R=/work/nvme/bivj/jma20/w121_mhd_0929
cd $R
while :; do
  if [ -n "$(squeue -u $USER -h -p ghx4-interactive -o %i)" ]; then sleep 60; continue; fi
  moved=0
  while read arm jid nr lim; do
    [ "$(squeue -h -j $jid -o %T 2>/dev/null)" = PENDING ] || continue
    o=""; [ $nr = 1 ] && o="--ntasks-per-node=1 --gpus-per-node=1 --mem=100G"
    new=$(sbatch --parsable -p ghx4-interactive -t $lim $o -J w$arm run_arm_deltaai.sub $R/runs/$arm 2 $nr 2>&1)
    if [[ $new =~ ^[0-9]+$ ]]; then
      scancel $jid
      echo "$(date '+%F %T') $arm: ghx4 $jid cancelled -> interactive $new ($lim)" | tee -a rotation.log
      sed -i "s/^$arm $jid /$arm $new /" rotation.txt; moved=1
    else
      echo "$(date '+%F %T') $arm: interactive submit failed: $new" | tee -a rotation.log
    fi
    break
  done < rotation.txt
  if [ $moved = 0 ] && [ -z "$(squeue -u $USER -h -p ghx4-interactive -o %i)" ]; then
    n=$(while read arm jid nr lim; do [ "$(squeue -h -j $jid -o %T 2>/dev/null)" = PENDING ] && echo x; done < rotation.txt | wc -l)
    [ $n = 0 ] && { echo "$(date '+%F %T') no pending ghx4 arm left, rotation done" | tee -a rotation.log; exit 0; }
  fi
  sleep 60
done
