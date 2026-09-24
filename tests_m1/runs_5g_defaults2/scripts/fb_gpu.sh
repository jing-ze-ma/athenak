#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/m1def2_0924/gpuf/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/m1def2_0924/gpuf/log.err.%j
#SBATCH -J m1def2f
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-defaults2: which h2fast lever raises the gas-Newton bracketed fallbacks under
# hesdirk2 (3-D He box 84x104x104, 4 blocks, 120 cycles, the runs_5f box), new binary
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/m1def2_0924; G=$W/gpuf; E=$W/bin/athena_new_box_gpu
md5sum $E
P="problem/vpert=1.0e-2"
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
T="rad_m1/implicit_tol=1.0e-11 rad_m1/implicit_lin_tol=1.0e-11 rad_m1/time2_lin_tol_fac=1"
declare -a A=(
 "E_D|" "E_L1|rad_m1/time2_lin_tol_fac=1" "E_S0|rad_m1/time2_one_pass_safety=0"
 "E_L1S0|rad_m1/time2_lin_tol_fac=1 rad_m1/time2_one_pass_safety=0"
 "E_OP0|rad_m1/implicit_one_pass=0" "E_T|$T"
 "V_D|$V" "V_L1|$V rad_m1/time2_lin_tol_fac=1" "V_S0|$V rad_m1/time2_one_pass_safety=0"
 "V_T|$V $T")
for a in "${A[@]}"; do
  IFS='|' read -r n ov <<< "$a"
  d=$G/$n; rm -rf $d; mkdir -p $d; cd $d
  echo "#### $n $(date +%T)"
  srun -n 1 $E -i $W/inp/box3d_fb.athinput -d . -t 00:02:30 time/nlim=120 $P $ov > run.log 2> run.err
  echo "rc=$? $(date +%T)"
done
python3 $W/scripts/fbcmp.py $G > $G/RESULTS.txt 2>&1
echo FB DONE
