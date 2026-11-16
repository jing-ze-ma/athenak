#!/bin/bash
# rt.sh NAME BIN DECK THREADS [overrides]: 16 cycles vs 8 + restart + 8, cmp of cycle-16 bins
N=$1; B=$2; DK=$3; T=$4; shift 4
D=$REG/rst/runs/$N; rm -rf $D; mkdir -p $D/f $D/h
export OMP_PROC_BIND=false OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=$T
(cd $D/f && $B -i $DK time/nlim=16 "$@" > a.log 2>&1)
(cd $D/h && $B -i $DK time/nlim=8 "$@" > a.log 2>&1 && r=$(ls rst/*.rst | tail -1) && $B -r $r time/nlim=16 > b.log 2>&1)
res=$(python3 $REG/rst/cmpbin.py $D/f $D/h 2>&1 | tail -1)
echo "RT $N:$res | $(grep -ih 'not bitwise\|warning: restart' $D/h/b.log | cut -c1-120 | tr '\n' ' ')"
