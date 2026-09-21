#!/bin/bash
# tests_mhd_rst/gate_hyd.sh -- the HYDRO general-EOS restart gate, for the no-regression
# check: the hydro half of this mechanism was closed by 6ae9ccbe and must stay closed.
# Same shape as gate.sh; see it and README.md.
set -u
BIN=$1
TAG=$2
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE" || exit 1
RUN="hyd_${TAG}"
rm -rf "${RUN}_A" "${RUN}_B"
mkdir -p "${RUN}_A" "${RUN}_B"
OVR="output3/dt=1.0e-30"
"$BIN" -i lw_hyd_geneos.athinput -d "${RUN}_A" time/nlim=20 $OVR > "${RUN}_A.log" 2>&1
"$BIN" -i lw_hyd_geneos.athinput -d "${RUN}_B" time/nlim=10 $OVR > "${RUN}_B1.log" 2>&1
RST=$(ls -1 "${RUN}_B"/rst/*.rst | tail -1)
"$BIN" -r "$RST" -d "${RUN}_B" time/nlim=20 > "${RUN}_B2.log" 2>&1
python3 cmpdumps.py "${RUN}_A/bin" "${RUN}_B/bin" | tee "${RUN}.cmp"
rm -rf "${RUN}_A/bin" "${RUN}_B/bin" "${RUN}_A/rst" "${RUN}_B/rst"
