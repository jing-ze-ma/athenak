#!/bin/bash
# wait for the CPU binaries, then gates_fuse.py + rst_mg.py
W=/viper/ptmp2/jinma/mgf_0925; T=/viper/ptmp2/jinma/wt_mgfuse/tests_m1
until grep -q athena_base_box_convection_cpu $W/b_basec.log 2>/dev/null; do sleep 60; done
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
python3 $T/runs_5n_mgfuse/gates_fuse.py run $W/bin/athena_base_box_convection_cpu $W/bin/athena_v1_box_convection_cpu $W/gates > $W/gates_run.log 2>&1
python3 $T/runs_5n_mgfuse/gates_fuse.py eval $W/gates > $W/gates_eval.txt 2>&1
ATHENAK_M1_DATA=/viper/ptmp2/jinma/faces_0924/m1data python3 $T/runs_5m_precond/rst_mg.py $W/bin/athena_v1_box_convection_cpu $W/rstmg > $W/rstmg.log 2>&1
echo CPU_GATES_DONE
