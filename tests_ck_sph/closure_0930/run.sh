#!/bin/bash
# usage: run.sh <bin> <name> <input> [key=val ...]   (login node, serial, nice)
BIN=$1; NAME=$2; IN=$(readlink -f $3); shift 3
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
D=/viper/ptmp2/jinma/cksph_test_0930/runs/$NAME
rm -rf $D; mkdir -p $D; cd $D
OMP_NUM_THREADS=1 nice -n 10 $BIN -i $IN -d $D time/nlim=2 time/tlim=1e30 \
  problem/ck_dump_file=$D/dump.txt \
  problem/ck_int_at_cut=false "$@" > log.txt 2>&1
echo "$NAME rc=$? $(grep -c FATAL log.txt) fatal"
