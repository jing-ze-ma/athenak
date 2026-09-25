#!/bin/bash -l
# usage: study.sh <list file> [parallel case/arm pairs]: st.py run for every line
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
cd $(dirname $(readlink -f $0))/..
xargs -P ${2:-8} -L 1 python3 st.py run < $1
echo STUDY DONE $1
