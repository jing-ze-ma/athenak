#!/bin/bash
# usage: evalall.sh <list file>: st.py eval for every (case, arm) line whose runs exist
cd /viper/ptmp2/jinma/wt_sporder2b/tests_m1/runs_5q_sporder2b
while read c a; do python3 st.py eval $c $a 2>&1 | grep -v "^Traceback\|^  File\|^    " ; done < $1
