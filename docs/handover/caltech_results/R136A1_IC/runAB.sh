#!/bin/bash -l
cd /resnick/groups/carnegie_poc/jingze/r136a1_1104/ic
export OPENBLAS_NUM_THREADS=1
module load gcc/13.2.0-gcc-13.2.0-w55nxkl hpcx/2.17.1/hpcx-ompi >/dev/null 2>&1
cd /home/jingze/ATHENAK/athenak/tools/ic_star
/resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python3 make_star_ic.py /resnick/groups/carnegie_poc/jingze/r136a1_1104/ic/r136a1B_noseam.json --out /resnick/groups/carnegie_poc/jingze/r136a1_1104/ic/AB
echo RC=$?
