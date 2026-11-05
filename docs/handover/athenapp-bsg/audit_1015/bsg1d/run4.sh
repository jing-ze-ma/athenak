#!/bin/bash
cd /raven/ptmp/jinma/athpp_audit_1015/bsg1d
(cd nb4 && ../../athena_nb/bin/athena -i athinput -d out > run.log 2>&1) &
(cd gv4 && timeout 600 ../../athena_bsg/bin/athena -i athinput -d out > run.log 2>&1) &
wait
