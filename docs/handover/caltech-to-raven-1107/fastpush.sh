#!/bin/bash
# fastpush.sh [K=3]  (Caltech login node, inside tmux, from $J; USER types the OTPs when prompted)
# Opens K extra INDEPENDENT connections to Raven (own gate hop each: gate OTP, then raven OTP), then copies
# RSG then HG as byte ranges spread over all connections (existing master + K new) into the .part names viper watches.
J=/resnick/groups/carnegie_poc/jingze; cd $J || exit 1
K=${1:-3}; R=/raven/ptmp/jinma
CP=(); m=$(ssh -G raven | awk '/^controlpath /{print $2}')
ssh -o ControlPath="$m" -O check raven 2>/dev/null && CP+=("$m")       # existing master (if alive)
for k in $(seq 2 $((K+1))); do
  p=~/.ssh/cmk$k-%r@%h:%p
  if ! ssh -o ControlPath=$p -O check raven 2>/dev/null; then
    echo ">>> connection $k: enter the GATE OTP, then the RAVEN OTP"
    ssh -o ControlMaster=yes -o ControlPersist=12h -o ControlPath=$p \
        -o ProxyCommand="ssh -o ControlMaster=no -o ControlPath=none -W %h:%p mpcdf-gate" -fN raven || { echo "connection $k failed"; continue; }
  fi
  CP+=("$p")
done
N=${#CP[@]}; echo "connections: $N"; [ $N -ge 1 ] || exit 1
S=$((N*4))                                                              # 4 streams per connection
pcp() { F=$1; D=$2; SZ=$(stat -c %s "$F"); MB=$(( SZ/1048576/S + 1 )); T0=$(date +%s)
  for i in $(seq 0 $((S-1))); do c=${CP[$((i % N))]}
    dd if="$F" bs=1M skip=$((i*MB)) count=$MB iflag=fullblock status=none \
      | ssh -o Compression=no -o ControlPath="$c" raven "dd of='$D' bs=1M seek=$((i*MB)) conv=notrunc iflag=fullblock status=none" &
  done; wait; dt=$(( $(date +%s)-T0 ))
  echo "$(date +%T) $(basename $F): $SZ B in ${dt}s = $(( SZ/1000000/(dt>0?dt:1) )) MB/s over $N connections x 4 streams; remote size $(ssh -o ControlPath="${CP[0]}" raven "stat -c %s '$D'")"; }
pcp rsg_E_1020/prodH/rst/rsg.01295.rst          $R/raven_gap_1105/rsgE/rst/rsg.01295.rst.part
pcp hgglob_1021/prodG/rst/hegiant.00139.rst      $R/raven_gap_1105/hg/rst/hegiant.00139.rst.part
echo "done; viper md5-checks the .part files and renames them"
