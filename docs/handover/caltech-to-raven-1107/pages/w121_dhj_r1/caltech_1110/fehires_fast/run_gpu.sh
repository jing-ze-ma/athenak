#!/bin/bash
# run_gpu.sh PLANET DUMP TAG OUTDIR [EOS] [extra key=val ...]
#   one trace for both cuts on the GPU -> OUTDIR/spec_<TAG>_pc6.npz, spec_<TAG>_pc0.npz, then ana.py on each
#   (+ the *_mx.json ana when planets/<PLANET>_mx.json exists, as run1.sh).  Set CUDA_VISIBLE_DEVICES first.
FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast; PY=/resnick/groups/carnegie_poc/jingze/fewind_1102/venv_fehgpu/bin/python
PJ=$1; D=$2; TAG=$3; OUT=$4; EOS=${5:-/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt}; shift 5 2>/dev/null || shift $#
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=8 NPROC=${NPROC:-8}
mkdir -p $OUT; cd $FF/code; S=$(date +%s)
$PY fe3d.py --planet planets/$PJ.json --eos $EOS --out $OUT $D tag=$TAG pcuts=1e-6,0 gpu=1 "$@" > $OUT/fe3d_$TAG.log 2>&1 || { echo FAIL fe3d $TAG; exit 1; }
echo "fe3d wall $(( $(date +%s)-S )) s" >> $OUT/fe3d_$TAG.log
# 11-03: the ana passes (pc6, pc0, and the _mx ones) run concurrently (CPU, ~10 s each); was serial
[ -f planets/${PJ}_mx.json ] && mkdir -p ${OUT}_mx
for pc in pc6 pc0; do
  ( $PY ana.py --planet planets/$PJ.json --out $OUT ${TAG}_$pc > $OUT/ana_${TAG}_$pc.log 2>&1 || echo FAIL ana ${TAG}_$pc ) &
  if [ -f planets/${PJ}_mx.json ]; then
    ln -sf $OUT/spec_${TAG}_$pc.npz ${OUT}_mx/spec_${TAG}_$pc.npz
    ( $PY ana.py --planet planets/${PJ}_mx.json --out ${OUT}_mx ${TAG}_$pc > $OUT/anamx_${TAG}_$pc.log 2>&1 || echo FAIL anamx ${TAG}_$pc ) &
  fi
done; wait; echo DONE $TAG
