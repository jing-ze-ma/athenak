#!/bin/bash
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103
srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c4 -J fe_mhd bash $O/worker3.sh 0 >> $O/logs/worker_g0.log 2>&1
srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c4 -J fe_bench bash -c "bash $O/dev/bench.sh 0 g1 FEHIRES_GATHER=1; bash $O/dev/bench.sh 0 g2 FEHIRES_GATHER=2" >> $O/logs/bench.log 2>&1
echo BENCHDONE >> $O/logs/bench.log
