#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
for a in Be0 Be1; do (cd $W/h2h/runs/$a && $W/athena_bc/bin/athena -i athinput -d outb > runb.log 2>&1) & done
wait
