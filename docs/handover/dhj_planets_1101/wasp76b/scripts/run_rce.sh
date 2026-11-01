#!/bin/bash
# usage: run_rce.sh <outdir> <G0> <AP> <T0.npz> [maxit]   (face-0 route-B RCE, rce_w121.py read-only)
PY=/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/bin/python3
source ${ENVF:-/viper/ptmp2/jinma/dhj_planets_1101/wasp76b/scripts/rce_env.sh}
mkdir -p $1 || exit 1
export RCE_G0=$2 RCE_AP=$3 RCE_OUT=$1/
cd $1 || exit 1
nice $PY /viper/ptmp2/jinma/wasp121_0925/rce_w121.py B 2 20 $4 ${5:-400} > $1/log.txt 2>&1
