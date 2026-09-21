#!/bin/bash -l
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
O=$R/tests_m1/runs_1cB
D=$O/t8a_beam_full; rm -rf $D; mkdir -p $D
(cd $D && $X -i $O/t8_beam.athinput > run.log 2>&1; echo "exit=$?" >> run.log)
D=$O/t8a_beam_rst; rm -rf $D; mkdir -p $D/rst
cp $O/t8a_beam_full/rst/m1_t8beam.00001.rst $D/rst/
(cd $D && $X -r rst/m1_t8beam.00001.rst > run.log 2>&1; echo "exit=$?" >> run.log)
echo T8A_BEAM_DONE
