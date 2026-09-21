#!/bin/bash
# one short dhj_ck_spherical run; $1 = dir, $2 = binary, rest = input overrides
BINDIR=$PWD
R=/viper/u2/jinma/ATHENAK/athenak
D=$1; BIN=$2; shift 2
rm -rf $D; mkdir -p $D
( cd $D && $BINDIR/$BIN -i $BINDIR/base.athinput \
    problem/ck_table=$R/data/exo_fms_ck/ck/Premixed_1x_g8_11.txt \
    problem/ck_data_dir=$R/data/exo_fms_ck \
    time/nlim=${NLIM:-20} output1/dt=1e-30 output3/dt=1e-30 output4/dt=1e30 \
    "$@" > run.log 2>&1 )
echo "$D exit=$?  $(grep 'zone-cycles' $D/run.log)"
