#!/bin/bash
# usage: bash SETUP.sh <files dir>   (absolute path; keep it for the whole test)
# gunzips the IC/tables into <files dir>, checks md5 (MD5SUMS = uncompressed files), writes the input
# with @FILES@ = <files dir>.  Can be run from anywhere.
set -e
H=$(cd "$(dirname "$0")" && pwd); F=$1
[ -n "$F" ] && [ "${F:0:1}" = / ] || { echo "usage: $0 <absolute files dir>"; exit 1; }
mkdir -p "$F"
for g in "$H"/*.txt.gz; do b=$(basename "$g" .gz); [ -f "$F/$b" ] || gunzip -c "$g" > "$F/$b"; done
(cd "$F" && md5sum -c "$H/MD5SUMS")
for i in "$H"/*.athinput.in; do sed "s#@FILES@#$F#g" "$i" > "$F/$(basename "$i" .in)"; done
grep -H "he_ic_file\|he_opac_table\|he_planck_table" "$F"/*.athinput
echo SETUP_OK
