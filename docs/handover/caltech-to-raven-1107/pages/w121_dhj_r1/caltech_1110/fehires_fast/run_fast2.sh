#!/bin/bash
# run_fast2.sh PLANET DUMP TAG OUTDIR [EOS] [extra key=val ...]   (11-03)
#   as run_gpu.sh but ONE cut (pcut 1e-6 bar -> OUTDIR/spec_<TAG>_pc6.npz, the observable) and prep=gpu
#   (ray geometry + nearest column + radial weights on the GPU, no CPU worker pool).  Set CUDA_VISIBLE_DEVICES.
FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast; PY=/resnick/groups/carnegie_poc/jingze/fewind_1102/venv_fehgpu/bin/python
PJ=$1; D=$2; TAG=$3; OUT=$4; EOS=${5:-/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt}; shift 5 2>/dev/null || shift $#
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $OUT; cd $FF/code; S=$(date +%s)
$PY fe3d.py --planet planets/$PJ.json --eos $EOS --out $OUT $D tag=$TAG pcuts=1e-6 tags=${TAG}_pc6 gpu=1 prep=gpu "$@" > $OUT/fe3d_$TAG.log 2>&1 || { echo FAIL fe3d $TAG; exit 1; }
echo "fe3d wall $(( $(date +%s)-S )) s" >> $OUT/fe3d_$TAG.log
$PY ana.py --planet planets/$PJ.json --out $OUT ${TAG}_pc6 > $OUT/ana_${TAG}_pc6.log 2>&1 &
if [ -f planets/${PJ}_mx.json ]; then
  mkdir -p ${OUT}_mx; ln -sf $OUT/spec_${TAG}_pc6.npz ${OUT}_mx/spec_${TAG}_pc6.npz
  $PY ana.py --planet planets/${PJ}_mx.json --out ${OUT}_mx ${TAG}_pc6 > $OUT/anamx_${TAG}_pc6.log 2>&1 || echo FAIL anamx ${TAG}_pc6
fi
wait; [ -s $OUT/ana_${TAG}_pc6.npz ] || echo FAIL ana ${TAG}_pc6
echo "total wall $(( $(date +%s)-S )) s" >> $OUT/fe3d_$TAG.log; echo DONE $TAG
