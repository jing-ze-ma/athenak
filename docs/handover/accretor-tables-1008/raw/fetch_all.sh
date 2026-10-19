#!/bin/bash
cd /viper/ptmp2/jinma/accretor_tables_1009 || exit
PT=/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/planck_tools
MIX="$(cat mixture_x0.70_z0.02.txt)"
f() { for a in 1 2 3; do date +%T; nice python3 $PT/fetch_tops.py "$MIX" "$@" && return 0; sleep 60; done; return 1; }
f 0.0005 1 1e-21 1e-14 36 raw/tops_lo_x0.70.dat > raw/fetch_lo.log 2>&1
f 0.0005 1 1e-14 1.0 71 raw/tops_gs98_x0.70_z0.02.dat > raw/fetch_hi.log 2>&1
f 1.0 10 1e-14 1.0 71 raw/tops_hiT_x0.70.dat > raw/fetch_hiT.log 2>&1
f 1.0 10 1e-21 1e-14 36 raw/tops_hiT_lo_x0.70.dat > raw/fetch_hiT_lo.log 2>&1
echo done > raw/fetch_all.done
