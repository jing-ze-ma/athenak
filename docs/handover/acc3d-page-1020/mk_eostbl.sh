#!/bin/bash
# Regenerate the run's general-EOS table dump for tgas.py / sthermo.py (only if you do not have a copy).
# The table is a pure function of the deck's <hydro> eos_* keys (eos = general, general_eos = table, X 0.70, Y 0.28, ionization on,
# no H2, no metal ionization, no radiation, logd -18..0, logT 3.5..7.6, dlog 0.02); AthenaK writes it when the EOS is built, i.e.
# at startup, so a 0-cycle run suffices. Expected: 901 x 1026 nodes, 25.7 MB, md5 e548d9c2777b8b63f979c024ab2dd564
# (Raven: /raven/ptmp/jinma/acc3d_1017/mk/eosdump/eos_x0.70_y0.28_logd-18.txt). A different md5 only matters if the eos_* keys differ.
# usage (inside a 1-GPU allocation):  BN=<accretor binary> DK=<substituted deck, e.g. your acc3d_prod_trc.athinput> \
#        DEST=<dir> LAUNCH="srun -n 1" ./mk_eostbl.sh
: ${BN:?} ${DK:?} ${DEST:?}
LAUNCH=${LAUNCH:-}
mkdir -p $DEST/eosdump_run && cd $DEST/eosdump_run
F=$DEST/eos_x0.70_y0.28_logd-18.txt
$LAUNCH $BN -i $DK hydro/eos_table_dump=$F time/nlim=0 problem/stream=false meshblock/nx3=512 -d $DEST/eosdump_run > run.log 2>&1
grep -m2 "General EOS" run.log
md5sum $F
