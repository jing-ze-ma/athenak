#!/bin/bash
# GPU cost breakdown of the production radiation scheme (DeltaAI port of delta-prof-1009, viper 10-09).
# usage (inside a 1-node allocation with 4 GH200, 4 ranks, 1 GPU and 1 MeshBlock per rank):
#   bash prof_cuda.sh <BIN = athena_hes_gpu_98835d99a16f> <AG Car files dir> <BSG files dir> <out dir> <KT .so>
# arms "name:case:mode:nlim"
#   case = agc  AG Car B production input (agcar_rcxB_ge_accel_st_mg_hr.athinput), FRESH (t = 0),
#               480x128x128 = 4 blocks 480x64x64, production outputs kept (hst 1e3 s, bins 1e4 s)
#          bsg  BSG half-range input (bsg_hr_dc5.athinput) on the reduced mesh of the mem batches
#               (256x128x128 = 4 blocks 256x64x64), bin/rst dumps switched off (the gate input dumps
#               every 5 cycles, production does not)
#   mode = plain  no tools: s/cycle
#          kt     KOKKOS_TOOLS_LIBS = simple-kernel-timer; *.dat per rank -> kp_reader *.txt
#          tmr    input copy with rad_m1/implicit_timers = 3 (the code's own fenced category timers)
# The kt arms run twice (short + long nlim); prof_group.py differences them so setup, IC, table
# reading and cycle 1 drop out (steady cycles only).
# env: ARMS, SRUN4, KPR (kp_reader; default next to the KT .so)
BIN=$1; FA=$2; FB=$3; S=$4; KT=$5
mkdir -p "$S" || exit 1
# DeltaAI binding (as agcar_dai_link.sh / mem_cuda*.sh): 72 Grace cores per rank, 1 GH200 per local rank
WRAP=$S/gpu_wrap.sh
printf '#!/bin/bash\nexport CUDA_VISIBLE_DEVICES=$SLURM_LOCALID\nexec "$@"\n' > "$WRAP"; chmod +x "$WRAP"
SRUN4=${SRUN4:-srun -n 4 -c 72 --cpu-bind=cores $WRAP}
KPR=${KPR:-$(dirname "$KT")/kp_reader}
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
IA=$FA/agcar_rcxB_ge_accel_st_mg_hr.athinput
IB=$FB/bsg_hr_dc5.athinput
GR="mesh/nx2=128 mesh/nx3=128 mesh/x2min=1.3089969389957472 mesh/x2max=1.832595714594046 mesh/x3max=0.5235987755982988 meshblock/nx2=64 meshblock/nx3=64"
KB="$GR output2/dcycle=1000000 output3/dcycle=1000000 output6/dcycle=1000000 output7/dcycle=1000000 output5/dt=1e30"
KA=""
[ -x "$BIN" ] || { echo "no binary $BIN"; exit 1; }
for f in "$IA" "$IB"; do [ -f "$f" ] || { echo "no input $f"; exit 1; }; done
[ -f "$KT" ] || { echo "no KT library $KT"; exit 1; }
[ -x "$KPR" ] || { echo "no kp_reader $KPR"; exit 1; }
echo "BIN $(md5sum $BIN)"; md5sum "$IA" "$IB"
go(){ # name case mode nlim
  D=$S/$1; mkdir -p $D && cd $D || exit 1
  if [ "$2" = agc ]; then I=$IA; K="$KA"; else I=$IB; K="$KB"; fi
  if [ "$3" = tmr ]; then
    T=$S/$(basename $I .athinput)_tmr.athinput
    [ -f $T ] || sed '/^<rad_m1>/a implicit_timers = 3' $I > $T
    grep -q '^implicit_timers = 3' $T || { echo "no timers line in $T"; return; }
    I=$T
  fi
  L=""; [ "$3" = kt ] && L=$KT
  echo "== $1 case=$2 mode=$3 nlim=$4 in=$(basename $I) keys=[$K] start $(date +%T)"; t0=$(date +%s)
  nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -l 1 > smi.log 2>&1 &
  SP=$!
  KOKKOS_TOOLS_LIBS=$L $SRUN4 $BIN -i $I -d $D $K time/nlim=$4 time/ndiag=1 \
    < /dev/null > out.log 2>&1
  rc=$?; kill $SP
  pk=$(awk -F, '{if($2>m[$1])m[$1]=$2} END{for(i in m) printf "gpu%s=%dMiB ", i, m[i]}' smi.log)
  # s/cycle from the elapsed= diag lines: mean over cycles 5..nlim, and the median cycle (no dumps)
  grep -o "elapsed=[0-9.e+-]* cycle=[0-9]*" out.log | sed 's/elapsed=//;s/ cycle=/ /' > cyc.txt
  spc=$(python3 - "$4" <<'EOF'
import sys
n = int(sys.argv[1]); t = {}
for l in open('cyc.txt'):
    a, b = l.split(); t[int(b)] = float(a)
k0 = 5 if n > 10 else 1
if n in t and k0 in t:
    d = sorted(t[c] - t[c - 1] for c in range(k0 + 1, n + 1) if c in t and c - 1 in t)
    print("mean(%d-%d)=%.4f median=%.4f" % (k0 + 1, n, (t[n] - t[k0]) / (n - k0), d[len(d) // 2]))
else:
    print("na")
EOF
)
  echo "rc=$rc fatal=$(grep -c FATAL out.log) oom=$(grep -ci 'out of memory\|failed to allocate' out.log)" \
       "nonconv=$(grep -c 'NON-CONVERGED after' out.log) nan=$(cat *.hst 2>/dev/null | grep -ci nan)" \
       "wall=$(( $(date +%s)-t0 ))s s/cycle: $spc peak: $pk $(date +%T)"
  grep -h "Picard iterations mean\|<rad_m1> timers" out.log | tail -3
  if [ "$3" = kt ]; then
    n=0; for f in *.dat; do [ -f "$f" ] || continue; $KPR $f > kt_${f%.dat}.txt 2>&1; n=$((n+1)); done
    echo "kt: $n rank files -> kt_*.txt"
  fi; }
for a in ${ARMS:-agc_plain:agc:plain:60 agc_kt10:agc:kt:10 agc_kt60:agc:kt:60 agc_tmr:agc:tmr:60 bsg_plain:bsg:plain:30 bsg_kt5:bsg:kt:5 bsg_kt30:bsg:kt:30 bsg_tmr:bsg:tmr:30}; do
  IFS=: read n c m l <<< "$a"
  go $n $c $m $l
  # smoke rule: the first arm is the smoke; stop on rc != 0 or FATAL
  if [ "$n" = agc_plain ] && { [ $rc -ne 0 ] || grep -q FATAL $S/$n/out.log; }; then
    echo "SMOKE FAILED ($n): stop"; exit 1
  fi
done
