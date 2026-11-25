#!/bin/bash
source /resnick/groups/carnegie_poc/jingze/e1s_1014/pagework/env.sh; source /resnick/groups/carnegie_poc/jingze/vis3d/venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/resnick/groups/carnegie_poc/jingze/vis3d/.mpl LP_NUM_THREADS=2
cd /resnick/groups/carnegie_poc/jingze/mov3d_1125; exec python w3d.py "$@"
