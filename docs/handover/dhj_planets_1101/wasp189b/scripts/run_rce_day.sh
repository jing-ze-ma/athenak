#!/bin/bash
# usage: run_rce_day.sh <outdir> <G0> <AP> <T0.npz> <TEQ_eff> [maxit]
# dayside 1-D estimate: the same global-mean RCE (rce_w121.py, mu0 sphere average) with the
# forcing raised: TEQ_eff = 2^(1/4) Teq (dayside-mean, no day-night redistribution) or
# sqrt(2) Teq (4x: substellar flux spread over the mu0 average).  Note rce_w121.py derives
# T_int from TEQ (Thorngren fit), so T_int differs from the global run (printed in log.txt).
PY=/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/bin/python3
source /viper/ptmp2/jinma/dhj_planets_1101/wasp189b/scripts/rce_env.sh
export RCE_TEQ=$5 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
mkdir -p $1 || exit 1
export RCE_G0=$2 RCE_AP=$3 RCE_OUT=$1/
cd $1 || exit 1
nice $PY /viper/ptmp2/jinma/wasp121_0925/rce_w121.py B 2 20 $4 ${6:-400} > $1/log.txt 2>&1
