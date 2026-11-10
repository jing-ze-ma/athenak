#!/bin/bash
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj/lock_hpc-sm-01-11_0
while kill -0 2284717 2>/dev/null; do sleep 5; done
srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c4 -J fe_mhd bash $O/worker3.sh 0 >> $O/logs/worker_g0.log 2>&1
rm -rf $L; echo "$(date +%T) GPU 0 released" >> $O/logs/hunt.log
