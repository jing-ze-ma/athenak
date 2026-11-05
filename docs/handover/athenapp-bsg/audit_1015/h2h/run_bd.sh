#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
for a in Bd0 Bd1; do (cd $W/h2h/runs/$a && $W/athena_bc/bin/athena -i athinput -d outb > runb.log 2>&1) & done
wait
