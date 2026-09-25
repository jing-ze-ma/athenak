#!/bin/bash
# usage: evalall.sh <list file>: st.py eval for every (case, arm) line whose runs exist
cd $(dirname $(readlink -f $0))/..
while read c a; do python3 st.py eval $c $a 2>&1 | grep -v "^Traceback\|^  File\|^    " ; done < $1
