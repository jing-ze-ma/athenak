export SRC=$HERE/src
export RUN=${RUN:-$R/prodF} BINDIR=${BINDIR:-$R/prodF/bin} EXTRA_RUN=${EXTRA_RUN:-$R/prodF} LOGDIR=${LOGDIR:-$R/logsF}
export OUT=${OUT:-$HERE/out} ICF=${ICF:-$R/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt}
export OLDRUN=""
export BASENAME=${BASENAME:-hegiant}
export RPH0=${RPH0:-7.94653e12}            # IC column header R_ph (ic_hegiant_14p5M_mk.txt)
export RIN=${RIN:-5.5656e11}               # deck <mesh>/x1min (inner wall), drawn dashed in the zooms
export TFLOOR=${TFLOOR:-2400}              # deck <hydro>/tfloor_kelvin
export JOBNAME=${JOBNAME:-hgGF} AGE=${AGE:-120}
export VISPY=${VISPY:-/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d}
export PY=${PY:-/resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python}   # numpy, scipy, matplotlib, h5py, PIL
export SITE=${SITE:-Caltech}
export LAYOUT=${LAYOUT:-"2 nodes × 4 H200 = 8 ranks, 24 MeshBlocks 704×64×64 (3 per GPU); 7 d links"}
export FFMPEG=${FFMPEG:-/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin/ffmpeg}
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MPLBACKEND=Agg
