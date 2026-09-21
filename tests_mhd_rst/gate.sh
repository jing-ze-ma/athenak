#!/bin/bash
# tests_mhd_rst/gate.sh -- the general-EOS MHD restart gate.
#
#   ./gate.sh <athena binary> <tag> [<meshblock nx1>] [<mpi ranks>]
#
# Runs the tabulated-EOS MHD fast magnetosonic wave (inputs/tests/
# linear_wave_mhd_geneos_table.athinput, built-in pgen "linear_wave") twice:
#   A: 20 cycles uninterrupted
#   B: 10 cycles, then restart from the cycle-10 restart file and run to cycle 20
# and compares the per-cycle binary dumps of the two.  A bitwise restart means every
# cycle of B equals the matching cycle of A byte for byte.
#
# The dumps are deleted at the end (inode quota); the comparison output is kept.
set -u
BIN=$1
TAG=$2
MBNX=${3:-64}
NRANK=${4:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE" || exit 1
RUN="run_${TAG}"
rm -rf "${RUN}_A" "${RUN}_B"
mkdir -p "${RUN}_A" "${RUN}_B"
LAUNCH=""
if [ "$NRANK" -gt 1 ]; then LAUNCH="mpirun -n $NRANK"; fi
OVR="meshblock/nx1=${MBNX} output3/dt=1.0e-30"
$LAUNCH "$BIN" -i lw_mhd_geneos.athinput -d "${RUN}_A" time/nlim=20 $OVR \
    > "${RUN}_A.log" 2>&1 || { echo "run A failed, see ${RUN}_A.log"; exit 1; }
$LAUNCH "$BIN" -i lw_mhd_geneos.athinput -d "${RUN}_B" time/nlim=10 $OVR \
    > "${RUN}_B1.log" 2>&1 || { echo "run B1 failed, see ${RUN}_B1.log"; exit 1; }
RST=$(ls -1 "${RUN}_B"/rst/*.rst | tail -1)
$LAUNCH "$BIN" -r "$RST" -d "${RUN}_B" time/nlim=20 \
    > "${RUN}_B2.log" 2>&1 || { echo "run B2 failed, see ${RUN}_B2.log"; exit 1; }
python3 cmpdumps.py "${RUN}_A/bin" "${RUN}_B/bin" | tee "${RUN}.cmp"
rm -rf "${RUN}_A/bin" "${RUN}_B/bin" "${RUN}_A/rst" "${RUN}_B/rst"
