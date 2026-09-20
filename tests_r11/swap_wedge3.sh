#!/bin/bash
# wait until wedge2 has passed 2.3 turnovers (or has ended), then take its nodes over
# with the entropy-seed arm wedge3
cd /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11
while squeue -h -j 11850256 | grep -q . ; do
  t=$(tail -1 wedge2/he4.hydro.hst | awk '{print ($1 > 10822) ? 1 : 0}')
  [ "$t" = "1" ] && break
  sleep 60
done
sbatch --parsable --export=ALL,OLDJOB=11850256 -o wedge3.log r11_wedge.sh wedge3 \
  time/cfl_number=0.15 problem/rt_top_vacuum=true problem/mlt_split_deposit=true \
  problem/cs_max=1.0e8 mesh/nx2=192 mesh/nx3=192 meshblock/nx2=24 meshblock/nx3=24 \
  problem/vpert_var=eint problem/vpert=3.0e-2 problem/vpert_sp_rand=true
