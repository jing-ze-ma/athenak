eval "$(grep -E '^(export|HERE=|R=)' /resnick/groups/carnegie_poc/jingze/hgglob_page/update.sh)"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
exec 9>/resnick/groups/carnegie_poc/jingze/hgglob_page/out/.update.lock; flock 9; cd /resnick/groups/carnegie_poc/jingze/hgglob_page/out
$PY $SRC/mk_movies.py data /resnick/groups/carnegie_poc/jingze/hgglob_page/out 2>&1 | grep -E "fps|encod|Error|rror" | tail -6
rsync -a /resnick/groups/carnegie_poc/jingze/hgglob_page/out/publish/movies/ /resnick/home/jingze/ATHENAK/pages/hggl_glob/movies/ --exclude "*_3d*"
