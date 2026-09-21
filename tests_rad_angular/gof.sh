#!/bin/bash
BINDIR=$PWD
R=/viper/u2/jinma/ATHENAK/athenak
D=$1; shift
rm -rf $D; mkdir -p $D
( cd $D && $BINDIR/athena_new -i $BINDIR/fatal.athinput \
    problem/ck_table=$R/data/exo_fms_ck/ck/Premixed_1x_g8_11.txt \
    problem/ck_data_dir=$R/data/exo_fms_ck "$@" time/nlim=0 > run.log 2>&1 )
echo "$D exit=$?"; grep -A4 "FATAL ERROR in.*conduction" $D/run.log | head -6
