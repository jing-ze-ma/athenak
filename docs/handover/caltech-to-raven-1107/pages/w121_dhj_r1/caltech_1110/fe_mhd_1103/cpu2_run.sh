#!/bin/bash
srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c8 --gres=none -J fe_mhd_cpu2 bash /resnick/groups/carnegie_poc/jingze/fe_mhd_1103/worker2.sh cpu stop_cpu2 >> /resnick/groups/carnegie_poc/jingze/fe_mhd_1103/logs/worker_cpu2.log 2>&1
