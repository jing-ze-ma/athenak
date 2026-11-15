#!/bin/bash
export OMP_NUM_THREADS=1 PYTHONPATH=/raven/ptmp/jinma/idort1020/page/py
cd /raven/ptmp/jinma/idort1020/page
/u/jinma/miniconda3/envs/magritte/bin/python figsrc/extract.py cb cbg shd hcyl nlte dif eq shock cloud marsh pipe g5 sp_core sp_nlte sp_shd sp_pinc sp_dif sp_relax sp_marsh sp_wave sp_patm
echo EXT_DONE
