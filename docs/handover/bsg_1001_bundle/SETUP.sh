#!/bin/bash
# usage: bash SETUP.sh <run dir>
# Writes runnable copies of the two BSG arm-2 inputs into <run dir>, with @BUNDLE@ replaced by the
# absolute path of this bundle directory (IC file + TOPS Rosseland/Planck tables live here).
set -e
B=$(cd "$(dirname "$0")" && pwd)
D=${1:?usage: SETUP.sh <run dir>}; mkdir -p "$D"
md5sum -c "$B/MD5SUMS"
for f in bsg3d_arm2.athinput bsg_col_arm2.athinput; do
  sed "s#@BUNDLE@#$B#g" "$B/$f" > "$D/$f"
done
grep -H 'he_ic_file\|he_opac_table\|he_planck_table' "$D"/*.athinput
