#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
(cd $W/h2h/runs/A_tc4 && $W/athena_bsg/bin/athena -i athinput -d out > run.log 2>&1) &
(cd $W/h2h/runs/A_nf4 && $W/athena_nf/bin/athena -i athinput -d out > run.log 2>&1) &
wait
