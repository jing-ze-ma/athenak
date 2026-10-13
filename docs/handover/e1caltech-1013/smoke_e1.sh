#!/bin/bash -l
#SBATCH -J agcE1smk
# e1caltech-1013 smoke_e1.sh SMOKE_DIR BIN BIN_MD5 DECK [LAUNCH]: PS's e1final recipe (Raven 31120638 / Viper 12172704-06, BITWISE)
# on 1 node x 2 H200, everything through the PRODUCTION link script link_e1.sh (same launch, same hook, same layout):
#   C    : fresh, 60 cycles                       (= the fresh-start smoke)
#   L    : link 1 fresh 30 cycles (ends with the link-end rst), then link 2 auto-restarts from the newest valid rst (+30 cycles)
#          (= restart BITWISE C 60 vs L 30+30, and the cadence check: every scheduled hst row / dump, none missing or extra)
#   Coff : C with rad_m1/vet_gd_halo_ipc=false (speed-only key: must be BITWISE vs C)
# Test-only keys (all exist in the deck): ndiag 1, hst 4e3, bins 1.5e4 (8 dumps per stem over 60 cycles), log output off.
# Submit: sbatch --export=ALL,E1C_DIR=<this dir> <your BSG flags for 1 node x 2 H200> --ntasks=2 -t 01:00:00 smoke_e1.sh ... Then: bash smoke_gate.sh SMOKE_DIR.
T=$1; BIN=$2; MD5=$3; DECK=$4; LAUNCH=${5:-}
W=${E1C_DIR:?export E1C_DIR=<dir of these scripts> (sbatch --export=ALL,E1C_DIR=...): under sbatch \$0 is a spooled copy}
[ -n "$DECK" ] || { echo "usage: smoke_e1.sh SMOKE_DIR BIN BIN_MD5 DECK [LAUNCH]"; exit 1; }
TK="time/ndiag=1 output1/dt=4e3 output2/dt=1.5e4 output3/dt=1.5e4 output4/dt=1e30"
mk() {  # dir ncyc extra-keys
  mkdir -p $T/$1
  cat > $T/$1/run.cfg <<CFG
E1C_DIR=$W
BIN=$BIN
BIN_MD5=$MD5
IN=$DECK
KEYS="$TK $3"
TLIM=2.4192e7
STOP_NC=1
E1C_NCYC=$2
${LAUNCH:+LAUNCH="$LAUNCH"}
CFG
}
mk C 60 ""; mk L 30 ""; mk Coff 60 "rad_m1/vet_gd_halo_ipc=false"
echo "== C";    bash -l $W/link_e1.sh $T/C
echo "== L1";   bash -l $W/link_e1.sh $T/L; cp -p $T/L/run.log.$SLURM_JOB_ID $T/L/run.log.L1; ls -l $T/L/rst
echo "== L2";   bash -l $W/link_e1.sh $T/L; ls -l $T/L/rst; cat $T/L/rst_prune.log 2>/dev/null
echo "== Coff"; bash -l $W/link_e1.sh $T/Coff
for d in C L Coff; do
  f=$T/$d/run.log.$SLURM_JOB_ID
  echo "--- $d: $(grep -c FATAL $f) FATAL; $(grep -m1 -E 'vet_gd_halo_ipc:' $f)"
  grep -E "implicit transport: .*Picard" $f | tail -1
  grep -E "he_ic_equil" $f | tail -2
  grep -E "^elapsed=" $f | sed -n "11p;\$p"   # s/cycle = (elapsed_last - elapsed_c10)/(cycles)
done
BINCONV_DIR=$W python3 $W/cmplast.py $T/C $T/Coff | tail -1
for h in $(cd $T/C; ls *.hst); do cmp -s <(grep -v '^#' $T/C/$h) <(grep -v '^#' $T/Coff/$h) && echo "IPC-OFF HST $h BITWISE" || echo "IPC-OFF HST $h DIFF"; done
nvidia-smi --query-gpu=index,memory.used --format=csv 2>/dev/null | head -3
echo E1_SMOKE_DONE
