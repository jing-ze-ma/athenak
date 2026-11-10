#!/bin/bash
srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c8 --gres=none -J fe_mhd_cpu bash /resnick/groups/carnegie_poc/jingze/fe_mhd_1103/worker.sh cpu >> /resnick/groups/carnegie_poc/jingze/fe_mhd_1103/logs/worker_cpu.log 2>&1
