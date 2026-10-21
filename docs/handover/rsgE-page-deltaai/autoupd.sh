#!/bin/bash
# autoupd.sh: rebuild the RSG E page every 45 min while the run is alive (niced, login node). Stop: touch STOP_autoupd. Log: autoupd.log
H=/work/nvme/bivj/jma20/rsgE_page
while [ ! -f $H/STOP_autoupd ]; do
  if squeue -h -u jma20 -n rsgE21 -t R 2>/dev/null | grep -q . || [ -n "$(find /projects/bivj/jma20/rsg_E_1021/run/bin -name '*.bin' -mmin -60 2>/dev/null | head -1)" ]; then
    flock -n $H/.lock nice -n 10 $H/update.sh > $H/autoupd.last 2>&1 && echo "$(date +%F_%T) ok" >> $H/autoupd.log || echo "$(date +%F_%T) FAIL (autoupd.last)" >> $H/autoupd.log
  fi
  sleep 2700
done
