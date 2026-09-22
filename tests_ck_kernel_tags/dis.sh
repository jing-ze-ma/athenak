#!/bin/bash -l
# Instruction / scratch / global-load counts per rt_chain_ck instantiation.
# Usage: ./dis.sh <deep_hot_jupiter_rt.cpp.o or athena binary>
# (same recipe as bench/bisect_cost/dis5.sh, with the tag-aware counter)
set -e
OBJ=$1
W=$(mktemp -d)
LLVM=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4/llvm/bin
$LLVM/llvm-objcopy --dump-section=.hip_fatbin=$W/fat.bin $OBJ /dev/null
$LLVM/clang-offload-bundler --type=o --unbundle --input=$W/fat.bin \
  --targets=hipv4-amdgcn-amd-amdhsa--gfx942:xnack+ --output=$W/co.o 2>/dev/null || \
$LLVM/clang-offload-bundler --type=o --unbundle --input=$W/fat.bin \
  --targets=hipv4-amdgcn-amd-amdhsa--gfx942 --output=$W/co.o
$LLVM/llvm-objdump -d $W/co.o > $W/d.txt
python3 "$(dirname "$0")/count.py" $W/d.txt
rm -rf $W
