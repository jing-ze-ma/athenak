#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
for a in L_tc4 L_Bd1 L_chk; do (cd $W/h2h/runs/$a && $W/athena_led/bin/athena -i athinput -d out > run.log 2>&1) & done
wait
