#!/bin/bash -l
# m1-h2acc: everything on the login node after the builds.  logs in $W/*.log
W=/viper/ptmp2/jinma/h2acc_0924; S=$W/scripts; cd $W
for b in base_none new_none exp_none base_box new_box; do while [ ! -x bin/athena_${b}_cpu ]; do sleep 60; done; done
( $S/gate_cpu.sh A B D > gate_ABD.log 2>&1 ) &
( $S/val.sh TV TV2 > val_TV.log 2>&1; $S/val.sh T1 > val_T1.log 2>&1 ) &
( $S/ss.sh > ss.log 2>&1 ) &
( for a in base new exp; do $S/shrun.sh ${a}_none sh_$a & done; $S/shrun.sh new_none sh_new_be rad_m1/time_scheme=be & wait; echo SH DONE > sh.done ) &
wait
( $S/rw.sh $W/lists/g1v.txt 16 new g1v_new > rw_g1v_new.log 2>&1; $S/rw.sh $W/lists/g1v.txt 16 base g1v_base > rw_g1v_base.log 2>&1 ) &
( $S/rw.sh $W/lists/g1x.txt 16 exp g1x > rw_g1x.log 2>&1 ) &
( $S/rw.sh $W/lists/stx.txt 12 exp stx > rw_stx.log 2>&1; $S/rw.sh $W/lists/st.txt 12 new st > rw_st.log 2>&1 ) &
( module purge; module load gcc/14 openmpi/5.0; python3 $S/runlist.py $W/lists/lmx.txt $W/le $W/bin/athena_exp_none_cpu 14 > le_x.log 2>&1; python3 $S/runlist.py $W/lists/len.txt $W/le $W/bin/athena_new_none_cpu 3 > le_n.log 2>&1 ) &
wait
echo ALL DONE > $W/all.done
