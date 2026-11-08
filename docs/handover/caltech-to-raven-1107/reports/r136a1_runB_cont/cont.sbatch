#!/bin/bash -l
#SBATCH -A carnegie_poc -p gpu -N1 -n2 -c 8 --gres=gpu:nvidia_h200:2 --mem=230G -t 7-00:00:00
#SBATCH --exclude=hpc-sm-01-09 -J r136B2 -o /resnick/groups/carnegie_poc/jingze/r136a1_1104/runB/cont/slurm.%j.out
#SBATCH --dependency=afterany:5544871
# R136a1 B continuation (USER 11-05 "make r136a1 longer"): restart runB/prod from its newest rst, same binary/deck/env/STOP rules as handover_run.sh prod loop.
J=/resnick/groups/carnegie_poc/jingze; R=$J/r136a1_1104; W=$R/runB; L=$W/launch.log
BIN=$J/builds_rstfixx/bin/athena_gpu_he_star_m1_prod_dc2dc20c_fmadoff; INFO=$J/bsg_1009/files/rst_info.py
# USER rule (since 11-06): name the poslim key on every restart (it is also carried in the rst params)
KEY="rad_m1/time2_adm_poslim=true"
[ -f $W/prod/STOP ] && { echo "$(date) CONT $SLURM_JOB_ID: prod/STOP exists, not continuing" >> $L; exit 1; }
pgrep -u $USER -f "athena_gpu_he_star_m1_prod_dc2dc20c_fmadoff.*$W/prod" > /dev/null && { echo "$(date) CONT: prod still running elsewhere, exit" >> $L; exit 1; }
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export OMPI_MCA_io=romio321 OMP_NUM_THREADS=1 HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1 UCX_TLS=^xpmem OMPI_MCA_coll_hcoll_enable=0
bad() { grep -aq "FATAL" $1 && return 0; grep -aqE "NON-CONVERGED=[1-9]" $1 && return 0; tail -n 60 $1 | grep -aqiw nan && return 0; return 1; }
cd $W/prod
for TRY in 0 1 2; do
  RS=$(ls rst/*.rst | sort | tail -1); INF=$(python3 $INFO $RS) || exit 1
  mv run.log run.$(date +%m%d_%H%M%S).log 2>/dev/null
  echo "$(date) CONT start (try $TRY) job $SLURM_JOB_ID on $(hostname) -r $RS -t 167:00:00 $KEY" >> $L
  srun -n 2 -c 8 --cpu-bind=none --mpi=pmix $J/r1planets_1101/tools/gpuwrap.sh 0 $BIN -r $RS $(echo "$INF" | tail -1) -d $W/prod -t 167:00:00 time/tlim=1.728e6 $KEY < /dev/null > run.log 2>&1; rc=$?
  echo "$(date) CONT rc $rc, newest rst $(ls rst/*.rst | sort | tail -1)" >> $L
  if bad run.log; then echo "$(date) CONT ended with FATAL/NON-CONVERGED/nan: STOP" | tee $W/prod/STOP >> $L; break; fi
  if [ $rc -eq 139 ] && [ $TRY -lt 2 ]; then sleep 30; continue; fi
  break
done
