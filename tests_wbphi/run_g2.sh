#!/bin/bash -l
# G2 of the wb-phieff branch: the three static-column arms, plus the restart check.
# Run from tests_wbphi/ with the box_convection binary as $1.
A=${1:?usage: run_g2.sh /path/to/athena}
for arm in A B C; do
  D=g2/$arm; rm -rf $D; mkdir -p $D
  (cd $D && $A -i ../../g2_arm${arm}.athinput > run.log 2>&1; echo "arm $arm exit=$?")
done
# G3: restart arm A from its mid-run file and compare the overlapping hst lines
rm -rf g2/Arst; mkdir -p g2/Arst; cp g2/A/rst/wb1d.00001.rst g2/Arst/
(cd g2/Arst && $A -r wb1d.00001.rst > run.log 2>&1; echo "arm A restart exit=$?")
python3 analyze_g2.py A B C
