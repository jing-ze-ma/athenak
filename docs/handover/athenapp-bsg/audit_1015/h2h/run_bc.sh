#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
R=$W/h2h/runs
(cd $R/Bchk && $W/athena_bc/bin/athena -i athinput -d out > run.log 2>&1) &
for a in Ba0 Ba1 Bb0 Bb1; do (cd $R/$a && $W/athena_bc/bin/athena -i athinput -d out > run.log 2>&1) & done
for a in Bc0 Bc1; do (cd $R/$a && $W/athena_bcnf/bin/athena -i athinput -d out > run.log 2>&1) & done
wait
