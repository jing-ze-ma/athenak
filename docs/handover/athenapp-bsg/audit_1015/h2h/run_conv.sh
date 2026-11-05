#!/bin/bash
W=/raven/ptmp/jinma/athpp_audit_1015
cd $W/h2h/runs/C_Bd1 && $W/athena_led/bin/athena -i athinput -d out > run.log 2>&1
