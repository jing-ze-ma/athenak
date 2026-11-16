#!/bin/bash
# rtn.sh NAME BIN DECK THREADS NTOT NRST: NTOT cycles straight vs NRST + restart + rest; cmp of the final bins
# (the deck writes its bins at cycle NTOT and a restart at NRST)
N=$1; B=$2; DK=$3; T=$4; NT=$5; NR=$6
D=$REG/rst/runs/$N; rm -rf $D; mkdir -p $D/f $D/h
export OMP_PROC_BIND=false OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=$T
(cd $D/f && $B -i $DK time/nlim=$NT > a.log 2>&1)
(cd $D/h && $B -i $DK time/nlim=$NR > a.log 2>&1 && r=$(ls rst/*.rst | tail -1) && $B -r $r time/nlim=$NT > b.log 2>&1)
res=$(python3 $REG/rst/cmpbin.py $D/f $D/h 2>&1 | tail -1)
echo "RTN $N:$res | ring: $(grep -ih 'rot_ring\|ring' $D/f/a.log | head -1 | cut -c1-100) | $(grep -ih 'not bitwise\|warning: restart' $D/h/b.log | cut -c1-120 | tr '\n' ' ')"
