#!/bin/bash -l
# T8a: a restart at mid-run must reproduce the uninterrupted run BITWISE.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
O=$R/tests_m1/runs_1cB
# --- T4 dynamic (coupled, moving, sub-cycled), 512 cells
D=$O/t8a_full; rm -rf $D; mkdir -p $D
(cd $D && $X -i $O/t8_advect.athinput > run.log 2>&1; echo "exit=$?" >> run.log)
D=$O/t8a_rst; rm -rf $D; mkdir -p $D/rst
cp $O/t8a_full/rst/m1_t8.00001.rst $D/rst/
(cd $D && $X -r rst/m1_t8.00001.rst > run.log 2>&1; echo "exit=$?" >> run.log)
# --- beam (radiation alone, no hydro), 128^2
D=$O/t8a_beam_full; rm -rf $D; mkdir -p $D
(cd $D && $X -i $O/t8_beam.athinput > run.log 2>&1; echo "exit=$?" >> run.log)
D=$O/t8a_beam_rst; rm -rf $D; mkdir -p $D/rst
cp $O/t8a_beam_full/rst/m1_t8beam.00001.rst $D/rst/
(cd $D && $X -r rst/m1_t8beam.00001.rst > run.log 2>&1; echo "exit=$?" >> run.log)
echo T8A_DONE
