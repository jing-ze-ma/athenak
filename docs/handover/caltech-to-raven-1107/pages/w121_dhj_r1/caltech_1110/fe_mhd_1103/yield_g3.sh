#!/bin/bash
# yield_g3.sh (coordinator 11-03): let the GPU 3 worker finish its current dump, stop it, release lock 3, wait 20 min, resume hunting.
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj/lock_hpc-sm-01-11_3
kill 981316                                   # hunter for GPU 3 (blocked in srun; its srun keeps running)
n0=$(grep -c " s [0-9]" $O/logs/worker_g3.log)
until [ $(grep -c " s [0-9]" $O/logs/worker_g3.log) -gt $n0 ]; do sleep 5; done
kill 981543 981544 2>/dev/null; sleep 15
run=$(timeout 30 ssh -n hpc-sm-01-11 "pgrep -u jingze -af fe3d.py" | grep -o 'tag=[^ ]*' | sed 's/tag=//')
for c in $(ls $O/claim); do
  [ -s $O/out/ana_${c}_pc0.npz ] && continue; echo "$run" | grep -qx "$c" && continue
  rm -f $O/out/*${c}*; rmdir $O/claim/$c; echo "released claim $c"
done
rm -rf $L; echo "$(date +%T) GPU 3 released for layout audit (20 min)" >> $O/logs/hunt.log
sleep 1200
echo "$(date +%T) resume hunting" >> $O/logs/hunt.log
exec bash $O/gpu_hunt3.sh
