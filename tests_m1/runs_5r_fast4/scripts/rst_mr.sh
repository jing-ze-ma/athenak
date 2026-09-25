#!/bin/bash
# MR restart gate: 14 cycles straight vs 7 + restart + 7 (k = $2), CPU, 4 ranks
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
EXE=$1; K=$2; W=/viper/ptmp2/jinma/fast4_0925; R=$W/cpu/rst_$(echo $K | tr " /=" "___"); rm -rf $R; mkdir -p $R/a $R/b
G="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=14 problem/vpert=1.0e-2 rad_m1/implicit_mr_every=$K"
cd $R/a; mpirun --oversubscribe -np 4 $EXE -i $W/inp/box_rst.athinput -d . $G > log.txt 2>&1
cd $R/b; mpirun --oversubscribe -np 4 $EXE -r $R/a/rst/m1slab.00001.rst -d . > log.txt 2>&1
ls $R/a/rst $R/b/rst
python3 - <<PY
a=open('$R/a/rst/m1slab.00002.rst','rb').read(); b=open('$R/b/rst/m1slab.00002.rst','rb').read()
ia=a.find(b'<par_end>'); ib=b.find(b'<par_end>')
print('K=$K bitwise after par_end:', a[ia:]==b[ib:], len(a)-ia, len(b)-ib)
PY
