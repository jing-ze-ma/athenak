#!/bin/bash -l
# vet_col atmosphere TRANSIENT (the T-S4 atmosphere, pure scattering rho kappa_s =
# 100 (r/r_in)^-2, from the Eddington start, c = 100): E after N = 400 steps of a FIXED
# dt = 2e-4 (t = 0.08, two light-crossing times), n = 32 .. 256 on the uniform grid.
# old = ref binary (surface_q), new = <tag> with implicit_marshak_face = linear and
# vet_col_order2 = true; ARMS="o2only mlin" runs inp/ts4_<arm>.athinput instead.
# Self-convergence: atmt.py.  usage: [ARMS=...] atmt.sh <tag> [ncore]
W=/viper/ptmp2/jinma/sporder2_0925; t=$1; nc=${2:-32}
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
N=400; DT=2.0e-4
SC="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=100.0"
for n in 32 64 128 256; do
  cfl=$(python3 -c "print(repr($DT*100.0*$n/4.0))")
  for a in ${ARMS:-old new}; do
    b=$W/bin/athena_${t}_none_cpu; [ $a = old ] && b=$W/bin/athena_ref_none_cpu
    d=$W/cpu/atmt/${a}${n}; rm -rf $d; mkdir -p $d
    (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $b -i $W/inp/ts4_$a.athinput -d $d \
      mesh/nx1=$n meshblock/nx1=$n time/nlim=$N output1/dcycle=$N $SC \
      rad_m1/implicit_cfl=$cfl rad_m1/vet_col_ncore=$nc rad_m1/implicit_tol=1.0e-12 \
      rad_m1/implicit_lin_tol=1.0e-12 problem/atm_seed=0.0 > log.txt 2>&1) &
  done
done
wait
echo ATMT DONE
