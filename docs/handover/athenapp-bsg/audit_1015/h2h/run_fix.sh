#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
for a in F_tc4 F_Bd1; do (cd $W/h2h/runs/$a && $W/athena_led/bin/athena -i athinput -d out > run.log 2>&1) & done
wait
