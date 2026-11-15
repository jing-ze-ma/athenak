#!/bin/bash
export PYTHONPATH=/raven/ptmp/jinma/idort1020/page/py
cd /raven/ptmp/jinma/idort1020/page; /u/jinma/miniconda3/envs/magritte/bin/python figsrc/plot.py "$@"; echo PLT_DONE
