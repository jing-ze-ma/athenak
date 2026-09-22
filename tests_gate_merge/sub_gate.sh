#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/athenak/tests_gate_merge/gate.log
#SBATCH -e /viper/u2/jinma/ATHENAK/athenak/tests_gate_merge/gate.err
#SBATCH -J rg_gate
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
# Bitwise gate of the he4-presn-global merge (e688efd4) against pre-merge rt-integration
# (8841786d): inputs/hydro/red_giant_cs.athinput, 100 cycles, one GPU, the same input file
# (the pre-merge one) through both binaries.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
A=/viper/u2/jinma/ATHENAK/athenak
G=$A/tests_gate_merge
I=/viper/u2/jinma/ATHENAK/bench/wt_ref/inputs/hydro/red_giant_cs.athinput
O="time/nlim=100 problem/ck_data_dir=$A/data/exo_fms_ck problem/ck_table=$A/data/exo_fms_ck/ck/Premixed_1x_g8_11.txt problem/opac_table=$A/data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt"
for v in ref new; do
  rm -rf $G/reg_$v; mkdir -p $G/reg_$v; cd $G/reg_$v
  srun -n 1 $G/athena_$v -i $I $O > log.txt 2>&1
done
cd $G
for f in rg.hydro.hst column.txt; do
  cmp reg_ref/$f reg_new/$f && echo "IDENTICAL $f" || echo "DIFFER $f"
done
md5sum reg_*/rg.hydro.hst reg_*/column.txt > regression_md5.txt
echo GATE_DONE
