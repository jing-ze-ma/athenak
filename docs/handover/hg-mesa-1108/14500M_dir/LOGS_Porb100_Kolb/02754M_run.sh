#!/bin/bash

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 23:00:00
##SBATCH --mem 8000
##SBATCH -p production 
#SBATCH -o 02754M.o
#SBATCH -e 02754M.e

export OMP_NUM_THREADS=8
cd /panfs/ds09/fuller/samanthawu/binary_rundir/r15140_runs/binary_test_new2/binary_test_EOScorr/14500M_dir
./clean
./mk
echo "SLURM JOB ID: $SLURM_JOB_ID"
cp inlist_to_end_core_he_burn_header inlist1
echo "    - inlist: inlist_to_end_core_he_burn_header " >> testhub.yml
./rn
cp inlist_to_end_core_c_burn_header inlist1
echo "    - inlist: inlist_to_end_core_c_burn_header " >> testhub.yml
./re
cp inlist_to_lgTmax_header inlist1
echo "    - inlist: inlist_to_lgTmax_header " >> testhub.yml
cp inlist_to_lgTmax_pt1 inlist_to_lgTmax
./re
cp inlist_to_lgTmax_pt2 inlist_to_lgTmax
./re 
cp inlist_to_cc_header inlist1
echo "    - inlist: inlist_to_cc_header " >> testhub.yml
./re
