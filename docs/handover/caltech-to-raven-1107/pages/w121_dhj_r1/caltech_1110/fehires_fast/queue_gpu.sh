#!/bin/bash
# queue_gpu.sh: runs fehires_fast/queue_gpu.txt (planet dump basetag) on the GPU (CUDA_VISIBLE_DEVICES set by caller)
# -> fehires/out/spec_<tag>_pc{6,0}.npz + ana (+ out_mx ana), logs in fehires/logs; ana runs in the background
F=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires; FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast
PY=/resnick/groups/carnegie_poc/jingze/fewind_1102/venv_fehgpu/bin/python; EOS=/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=8 NPROC=${NPROC:-8}
cd $FF/code
while read PJ D T; do
  [ -s $F/out/spec_${T}_pc6.npz ] && [ -s $F/out/spec_${T}_pc0.npz ] && { echo "skip $T"; continue; }
  S=$(date +%s); echo "$(date) start $T gpu $CUDA_VISIBLE_DEVICES $(hostname)"
  $PY fe3d.py --planet planets/$PJ.json --eos $EOS --out $F/out $D tag=$T pcuts=1e-6,0 gpu=1 > $F/logs/fe3d_${T}_gpu.log 2>&1 || { echo "FAIL fe3d $T"; continue; }
  echo "fe3d wall $(( $(date +%s)-S )) s" >> $F/logs/fe3d_${T}_gpu.log; echo "$(date) fe3d done $T $(( $(date +%s)-S )) s"
  ( for pc in pc6 pc0; do
      $PY ana.py --planet planets/$PJ.json --out $F/out ${T}_$pc > $F/logs/ana_${T}_$pc.log 2>&1 || echo "FAIL ana ${T}_$pc"
      if [ -f planets/${PJ}_mx.json ]; then mkdir -p $F/out_mx; ln -sf $F/out/spec_${T}_$pc.npz $F/out_mx/spec_${T}_$pc.npz
        $PY ana.py --planet planets/${PJ}_mx.json --out $F/out_mx ${T}_$pc > $F/logs/anamx_${T}_$pc.log 2>&1 || echo "FAIL anamx ${T}_$pc"; fi
    done ) &
done < $FF/queue_gpu.txt
wait; echo "$(date) QUEUEDONE"
