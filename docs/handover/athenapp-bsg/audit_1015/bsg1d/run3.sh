#!/bin/bash
cd /raven/ptmp/jinma/athpp_audit_1015/bsg1d
for d in nf4 nf1; do (cd $d && ../../athena_nf/bin/athena -i athinput -d out > run.log 2>&1) & done; wait
