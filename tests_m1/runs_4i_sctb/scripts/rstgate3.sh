#!/bin/bash -l
# restart gate: tb3_Fx2r4 restarted from rst 00001 (t = 1) must reproduce the t = 2, 3
# outputs of the continuous run byte for byte, and the hst tail
cd /viper/ptmp2/jinma/sctb_0923/cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
d=rs_tb3_Fx2r4; rm -rf $d; mkdir -p $d; cd $d
nice -n 10 mpirun -np 4 --oversubscribe --bind-to none ../../bin/athena_tb3_cpu -r ../tb3_Fx2r4/rst/m1slab.00001.rst > log.txt 2>&1
echo "exit $?" >> log.txt; cd ..
for f in bin/m1slab.m1.00002.bin bin/m1slab.m1.00003.bin bin/m1slab.hydro_w.00003.bin rst/m1slab.00002.rst rst/m1slab.00003.rst; do
  cmp -s tb3_Fx2r4/$f $d/$f && echo "same $f" || echo "DIFF $f"; done
for h in m1slab.user.hst m1slab.hydro.hst; do
  n=$(grep -vc '^#' $d/$h); diff <(grep -v '^#' tb3_Fx2r4/$h | tail -$n) <(grep -v '^#' $d/$h) > /dev/null && echo "hst tail ($n lines) same $h" || echo "hst DIFF $h"; done
tail -1 $d/log.txt
