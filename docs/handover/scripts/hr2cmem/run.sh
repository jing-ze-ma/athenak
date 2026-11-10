#!/bin/bash
# run.sh CASE ARM MODE "levels" : make the decks (cart.py, plane-parallel CART_1D=1) and run them one after the other.
# env: BIN (athena binary), CART_H (run tree, default ./cart_runs), MPI (default "mpirun -np 1")
HERE=$(cd $(dirname $0) && pwd); C=$1; A=$2; M=$3; LV=$4
export CART_1D=${CART_1D:-1}; MPI=${MPI:-mpirun -np 1}
for l in $LV; do inp=$(python3 $HERE/cart.py deck $C "$A" $M $l) || exit 1; d=$(dirname $inp)
  OMP_NUM_THREADS=1 $MPI $BIN -i $inp -d $d > $d/log 2>&1; echo "$C $A $M $l rc=$?" | tee -a ${CART_H:-./cart_runs}/verdicts.txt
done
python3 $HERE/cart.py eval $C "$A" $M $LV
