#!/bin/bash -l
#SBATCH -J bsgLIsmk
# bsgli-caltech-1014 smoke_bsgli.sh SMOKE_DIR BIN BIN_MD5 DECK [LAUNCH]: the production smoke + gate of the Caltech Linhao-IC BSG copy,
# 1 node x 2 H200, everything through the PRODUCTION link script link_bsgli.sh (same launch, same hook, same layout 2 ranks x 8 blocks):
#   C    : fresh, 60 cycles                       (= the fresh-start smoke)
#   L    : link 1 fresh 30 cycles (ends with the link-end rst), then link 2 auto-restarts from the newest valid rst (+30 cycles)
#          (= restart BITWISE C 60 vs L 30+30, and the cadence check: every scheduled hst row / dump, none missing or extra)
#   Coff : C with rad_m1/vet_gd_halo_ipc=false (speed-only key: must be BITWISE vs C, else halo_ipc is not kept)
# Then it runs smoke_gate.sh SMOKE_DIR ON THE COMPUTE NODE: GATE PASS writes SMOKE_DIR/VERIFIED and the job exits 0 (the production
# link 1 is submitted with --dependency=afterok on this job); GATE FAIL -> exit 1 (link 1 never starts).
# Test-only keys (all exist in the deck): ndiag 1, hst every 200 s, hydro_w/m1 bins every 1000 s, m1_vet/m1_face every 2500 s,
# log output off; rst as the deck (t = 0 and the link-end rst).
# Submit: sbatch --export=ALL,BLC_DIR=<this dir> <your BSG flags for 1 node x 2 H200> --ntasks=2 -t 01:15:00 -J bsgLIsmk smoke_bsgli.sh ...
T=$1; BIN=$2; MD5=$3; DECK=$4; LAUNCH=${5:-}
W=${BLC_DIR:?export BLC_DIR=<dir of these scripts> (sbatch --export=ALL,BLC_DIR=...): under sbatch \$0 is a spooled copy}
[ -n "$DECK" ] || { echo "usage: smoke_bsgli.sh SMOKE_DIR BIN BIN_MD5 DECK [LAUNCH]"; exit 1; }
TK="time/ndiag=1 output1/dt=200 output2/dt=1000 output3/dt=1000 output4/dt=1e30 output6/dt=2500 output7/dt=2500"
mk() {  # dir ncyc extra-keys   (LAUNCH written with printf %q: a multi-word LAUNCH keeps its quoting -- your E1 smoke 4452278 issue)
  mkdir -p $T/$1
  { echo "BLC_DIR=$W"; echo "BIN=$BIN"; echo "BIN_MD5=$MD5"; echo "IN=$DECK"; echo "KEYS=\"$TK $3\""
    echo "TLIM=4.96e6"; echo "STOP_NC=1"; echo "BLC_NCYC=$2"; [ -n "$LAUNCH" ] && printf 'LAUNCH=%q\n' "$LAUNCH"; } > $T/$1/run.cfg
}
mk C 60 ""; mk L 30 ""; mk Coff 60 "rad_m1/vet_gd_halo_ipc=false"
nvidia-smi --query-gpu=timestamp,index,memory.used --format=csv,noheader -l 5 > $T/gpumem.csv 2>/dev/null & MON=$!
echo "== C";    bash -l $W/link_bsgli.sh $T/C
echo "== L1";   bash -l $W/link_bsgli.sh $T/L; cp -p $T/L/run.log.$SLURM_JOB_ID $T/L/run.log.L1; ls -l $T/L/rst
echo "== L2";   bash -l $W/link_bsgli.sh $T/L; ls -l $T/L/rst; cat $T/L/rst_prune.log 2>/dev/null
echo "== Coff"; bash -l $W/link_bsgli.sh $T/Coff
kill $MON 2>/dev/null
for d in C L Coff; do
  f=$T/$d/run.log.$SLURM_JOB_ID
  echo "--- $d: $(grep -c FATAL $f) FATAL; $(grep -m1 -E 'vet_gd_halo_ipc:' $f)"
  grep -E "implicit transport: .*Picard" $f | tail -1
  grep -m2 -E "he_ic_linhao|he_bc_linhao" $f
  grep -E "cycle=" $f | awk '{for(i=1;i<=NF;i++) if($i ~ /^dt=/) print $i}' | sort | uniq -c | head -5
  grep -E "^elapsed=" $f | sed -n "11p;\$p"   # s/cycle = (elapsed_last - elapsed_c10)/(cycles)
done
echo "GPU peak MiB per GPU: $(awk -F', ' '{v=$3+0; if(v>m[$2]) m[$2]=v} END{for(i in m) printf "%s:%d ", i, m[i]}' $T/gpumem.csv 2>/dev/null)"
export BINCONV_DIR=$W
python3 $W/cmplast.py $T/C $T/Coff | tail -1
for h in $(cd $T/C; ls *.hst); do cmp -s <(grep -v '^#' $T/C/$h) <(grep -v '^#' $T/Coff/$h) && echo "IPC-OFF HST $h BITWISE" || echo "IPC-OFF HST $h DIFF"; done
echo "== GATE"; bash $W/smoke_gate.sh $T; g=$?
echo BSGLI_SMOKE_DONE gate_rc=$g
exit $g
