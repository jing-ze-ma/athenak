#!/bin/bash
# usage: bash SETUP.sh <dir holding these files>   -> checks MD5SUMS, writes agcar_shake{A,B}_ge.athinput with paths filled in
D=$(cd "$1" && pwd) || exit 1
cd "$D" || exit 1
md5sum -c MD5SUMS || { echo "MD5 FAIL"; exit 1; }
for c in A B; do sed "s#@AGCAR_DIR@#$D#g" agcar_shake${c}_ge.athinput.in > agcar_shake${c}_ge.athinput; done
grep -n "$D" agcar_shake?_ge.athinput
