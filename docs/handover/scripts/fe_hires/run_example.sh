#!/bin/bash
# WASP-121b reproduction on viper (10 G MHD s10_e15, dump 00166, rot 318).
# Reference (10-07): KpVsys full peak -4.26, segments mor -3.43 / eve -6.78 km/s.
D=$(cd "$(dirname "$0")" && pwd) || exit 1
OUT=${OUT:-$PWD/out_w121}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NPROC=${NPROC:-16}
python "$D/fe3d.py" --planet "$D/planets/wasp121b.json" \
  --eos /viper/ptmp2/jinma/w121_obs_1003/diag/3x/eos_table.txt --out "$OUT" \
  /viper/ptmp2/jinma/w121_obs_1003/mhd_final/data/s10_e15/dhj.mhd_w_bcc.00166.bin tag=s10_00166
python "$D/ana.py" --planet "$D/planets/wasp121b.json" --out "$OUT" s10_00166
