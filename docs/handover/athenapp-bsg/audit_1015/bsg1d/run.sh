#!/bin/bash
cd /raven/ptmp/jinma/athpp_audit_1015/bsg1d
for d in tc4 tc40 tc400 tc4n; do (cd $d && ../../athena_bsg/bin/athena -i athinput -d out > run.log 2>&1) & done; wait
