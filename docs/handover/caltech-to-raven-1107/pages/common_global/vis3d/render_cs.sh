#!/bin/bash
# CPU render job for a GLOBAL cubed-sphere run (cs3d.py): reduce new dumps (cached), render frames, encode, copy to publish/movies.
# usage: sbatch -A carnegie_poc -p expansion -c 8 --mem=64G -t 02:00:00 -J vis3d_RUN -o vis3d/RUN/job.%j.out render_cs.sh RUN [VIEW ...] [-- extra args]
V=/resnick/groups/carnegie_poc/jingze/vis3d
RUN=$1; shift
VIEWS=${VIEWS:-single layer}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=$V/.mpl PYTHONDONTWRITEBYTECODE=1
for v in $VIEWS; do
  /resnick/groups/carnegie_poc/jingze/venv_vtk/pvpy.sh $V/cs3d.py $RUN --view $v --workers ${WORKERS:-6} "$@"
done
