#!/bin/bash
# tests_cs_rst/gate_cs.sh -- the general-EOS CUBED-SPHERE restart gate.
#
#   ./gate_cs.sh <athena binary> <input file> <tag> [<meshblock nx2=nx3>] [<mpi ranks>]
#
# Same shape as tests_mhd_rst/gate.sh (and it reuses that directory's cmpdumps.py):
#   A: 20 cycles uninterrupted
#   B: 10 cycles, restart from the cycle-10 restart file, run to cycle 20
# and compares the per-cycle binary dumps bit for bit.  Dumps and restart files are
# deleted at the end (inode quota); the comparison table is kept.
set -u
BIN=$1
INP=$2
TAG=$3
MBNX=${4:-16}
NRANK=${5:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE" || exit 1
RUN="run_${TAG}"
rm -rf "${RUN}_A" "${RUN}_B"
mkdir -p "${RUN}_A" "${RUN}_B"
LAUNCH=""
if [ "$NRANK" -gt 1 ]; then LAUNCH="mpirun -n $NRANK"; fi
OVR="meshblock/nx2=${MBNX} meshblock/nx3=${MBNX} output3/dt=1.0e-30"
$LAUNCH "$BIN" -i "$INP" -d "${RUN}_A" time/nlim=20 $OVR \
    > "${RUN}_A.log" 2>&1 || { echo "run A failed, see ${RUN}_A.log"; exit 1; }
$LAUNCH "$BIN" -i "$INP" -d "${RUN}_B" time/nlim=10 $OVR \
    > "${RUN}_B1.log" 2>&1 || { echo "run B1 failed, see ${RUN}_B1.log"; exit 1; }
RST=$(ls -1 "${RUN}_B"/rst/*.rst | tail -1)
$LAUNCH "$BIN" -r "$RST" -d "${RUN}_B" time/nlim=20 \
    > "${RUN}_B2.log" 2>&1 || { echo "run B2 failed, see ${RUN}_B2.log"; exit 1; }
python3 ../tests_mhd_rst/cmpdumps.py "${RUN}_A/bin" "${RUN}_B/bin" | tee "${RUN}.cmp"
rm -rf "${RUN}_A/bin" "${RUN}_B/bin" "${RUN}_A/rst" "${RUN}_B/rst"
