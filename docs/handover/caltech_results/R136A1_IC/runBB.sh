#!/bin/bash -l
cd /resnick/groups/carnegie_poc/jingze/r136a1_1104/ic
export OPENBLAS_NUM_THREADS=1
/resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python3 /resnick/groups/carnegie_poc/jingze/r136a1_1104/tools/seam_only_B.py
echo RC=$?
