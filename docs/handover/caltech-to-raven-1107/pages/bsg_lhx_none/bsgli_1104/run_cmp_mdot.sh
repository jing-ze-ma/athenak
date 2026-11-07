#!/bin/bash -l
# run_cmp_mdot.sh: cmp_mdot.py inside alloc 5719508 (CPU only, 4 cores, 48 GB); output cmp_mdot.txt / cmp_mdot.png
W=/resnick/groups/carnegie_poc/jingze/bsgli_1104
srun --overlap --jobid=5719508 -w hpc-sm-01-05 -N1 -n1 -c 4 --mem=48G nice -n 10 bash -lc "source /etc/profile >/dev/null 2>&1; module purge; module load python/3.11.6-gcc-11.3.1-rphh4kv py-scipy/1.11.1-gcc-11.3.1-swvnmdi py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1; python3 $W/cmp_mdot.py $W/cmp_mdot.png" > $W/cmp_mdot.txt 2>&1
