#!/bin/bash
# cmds.txt: production R=1e6 on the ext2 node spacing (log T 0.025, log rho 0.05), then
# convergence R=1e4..3e5 and sensitivity runs at R=1e5 on the coarse grid (0.05 x 0.25)
P=/orion/ptmp/jinma/rsg_wind_1008/grains/venv/bin/python
S=/orion/ptmp/jinma/agcar_opac_1009/scripts/run_T.py
L=/orion/ptmp/jinma/agcar_opac_1009/work/logs
: > cmds.txt
for t in $(seq -f "%.3f" 3.45 0.025 4.301); do t2=$(printf "%.2f" $t); [ "$t" == "${t2}0" ] || t2=$t
  echo "LINES=mrg DLR=0.05 $P $S $t 1e6 _prod > $L/T${t}_prod.log 2>&1" >> cmds.txt; done
TS=$(seq -f "%.2f" 3.45 0.05 4.301)
for R in 1e4 3e4 1e5 3e5; do for t in $TS; do echo "LINES=mrg $P $S $t $R > $L/T${t}_$R.log 2>&1" >> cmds.txt; done; done
for t in $TS; do
  echo "LINES=gf08 $P $S $t 1e5 _gf08 > $L/T${t}_gf08.log 2>&1" >> cmds.txt
  echo "LINES=mrg WCAP=0.003 $P $S $t 1e5 _wcap003 > $L/T${t}_wcap003.log 2>&1" >> cmds.txt
  echo "LINES=mrg WCAP=0.03 $P $S $t 1e5 _wcap03 > $L/T${t}_wcap03.log 2>&1" >> cmds.txt
  echo "LINES=mrg XI_TURB=0 $P $S $t 1e5 _xi0 > $L/T${t}_xi0.log 2>&1" >> cmds.txt
  echo "LINES=mrg NMAX_H=15 $P $S $t 1e5 _nh15 > $L/T${t}_nh15.log 2>&1" >> cmds.txt
  echo "LINES=mrg NMAX_H=60 $P $S $t 1e5 _nh60 > $L/T${t}_nh60.log 2>&1" >> cmds.txt
  echo "LINES=mrg EPS=1e-12 $P $S $t 1e5 _eps12 > $L/T${t}_eps12.log 2>&1" >> cmds.txt
done
wc -l cmds.txt
