#!/bin/bash
# final CPU gates of m8: gates.py arms (box, slab; loose + tight) b0 vs m8 bitwise; restart
# bitwise; the runs_5q sp space-time set (f4 = m4 vs f8 = m8)
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/fast4_0925; G=/viper/ptmp2/jinma/wt_fast4/tests_m1/gates
cd $G
for t in b0 m8; do rm -rf $W/cpu/g_$t; python3 gates.py run $W/bin/athena_${t}_box_convection_cpu $W/cpu/g_$t box_a1 box_m2 slab_a1 slab_m2 > $W/cpu/g_$t.log 2>&1; done
for a in box_a1 box_m2 slab_a1 slab_m2; do for l in loose tight; do
  echo "== $a $l b0 vs m8: $(python3 cmp.py $W/cpu/g_b0/${a}_$l $W/cpu/g_m8/${a}_$l | tail -1)"; done; done
bash $W/scripts/rst_mr.sh $W/bin/athena_m8_box_convection_cpu 1 2>&1 | tail -1
cd $W/st; ln -sf $W/bin/athena_m8_none_cpu bin/athena_f8_none_cpu
grep -q "'f8'" st.py || sed -i "s/    'f7': ('f7', {}),/    'f7': ('f7', {}),\n    'f8': ('f8', {}),/" st.py
printf "pulse_u_T f8\nrsw_u_T f8\nmvr_u_T f8\npulse_u f8\natm_u f8\n" > list_f8.txt
module purge; module load gcc/14 openmpi/5.0
xargs -P 5 -L 1 nice python3 st.py run < list_f8.txt > run_f8.log 2>&1
for c in pulse_u_T rsw_u_T mvr_u_T pulse_u atm_u; do for a in f4 f8; do python3 st.py eval $c $a; done; done > RESULTS_f8.txt 2>&1
echo CPU8 DONE
