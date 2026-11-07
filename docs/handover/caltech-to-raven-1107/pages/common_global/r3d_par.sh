#!/bin/bash
# r3d_par.sh rsgE|hggl (10-23): full 3-D refresh in an expansion job: reductions RED_CAP 12, both views 24 workers, encode, stage.
# pagefast-1102 (old: r3d_par.sh.pre_fast_1102): RED_CAP 16 (one round for up to 16 new dumps on the 32-core step), and the
# layer view starts as soon as the single view has finished the (shared, cached) reductions instead of after its encode.
J=/resnick/groups/carnegie_poc/jingze; N=$1
case $N in rsgE) V3=rsgE_glob; PG=$J/rsgE_page;; hggl) V3=hggl_glob; PG=$J/hgglob_page;; esac
V=$J/vis3d; export MPLCONFIGDIR=$V/.mpl RED_CAP=${RED_CAP:-16} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; cd $V
L1=$V/$V3/job_render_single.log; L2=$V/$V3/job_render_layer.log
$J/venv_vtk/pvpy.sh $V/cs3d.py $V3 --view single --workers 24 > $L1 2>&1 & P1=$!
while kill -0 $P1 2>/dev/null && ! grep -q 'cached frames' $L1 2>/dev/null; do sleep 2; done
$J/venv_vtk/pvpy.sh $V/cs3d.py $V3 --view layer --workers 24 > $L2 2>&1 & P2=$!
wait $P1; r1=$?; wait $P2; r2=$?
echo "$N single rc=$r1 $(grep encoded $L1 | tail -1)"; echo "$N layer rc=$r2 $(grep encoded $L2 | tail -1)"
rsync -a $PG/out/publish/movies/ /resnick/home/jingze/ATHENAK/pages/${N}_glob/movies/ && ls -la /resnick/home/jingze/ATHENAK/pages/${N}_glob/movies/ | grep _3d
