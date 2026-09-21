#!/bin/bash -l
# G-recon (1c-B): the higher-order reconstructions.  pulse1d convergence, the beam, and
# the T3 AP rates.  Everything runs at <mesh>/nghost = 3 so that the four methods are
# compared on the same grid layout (plm at nghost 2 and 3 is bit-identical).
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
O=$R/tests_m1/runs_1cB
IP=$R/inputs/tests/rad_m1_pulse1d.athinput
IB=$R/inputs/tests/rad_m1_beam.athinput
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
run () { local n=$1; local IN=$2; shift 2; local d=$O/$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }

for rec in plm ppm4 ppmx wenoz; do
  for N in 32 64 128 256 512; do
    run p1d_${rec}_$N $IP mesh/nx1=$N meshblock/nx1=$N rad_m1/reconstruct=$rec &
  done
  run beam_$rec $IB mesh/nghost=3 rad_m1/reconstruct=$rec &
done
# plm at nghost 2, the 1c-A layout, to confirm nghost does not change the answer
run p1d_plm_ng2_128 $IP mesh/nghost=2 mesh/nx1=128 meshblock/nx1=128 &
wait
# T3: the AP rates for the two new high-order methods (plm and dc are in runs_1c)
for rec in ppm4 ppmx wenoz plm; do
  run t3_${rec}_tau1e1 $IT mesh/nghost=3 rad_m1/kappa_s=1280 rad_m1/reconstruct=$rec \
      time/tlim=0.5 output1/dt=0.05 &
  run t3_${rec}_tau1e3 $IT mesh/nghost=3 rad_m1/kappa_s=128000 rad_m1/reconstruct=$rec \
      time/tlim=12.0 output1/dt=1.2 &
  run t3_${rec}_tau1e6 $IT mesh/nghost=3 rad_m1/kappa_s=128000000 \
      rad_m1/reconstruct=$rec time/tlim=3000.0 output1/dt=300.0 &
  run t3_${rec}_nyq $IT mesh/nghost=3 rad_m1/kappa_s=128000 rad_m1/reconstruct=$rec \
      problem/nyquist_amp=1.0e-3 time/tlim=20.0 output1/dt=2.0 &
  for N in 64 128 256 512; do
    K=$(python3 -c "print(1280.0*$N/128)")
    run t3n_${rec}_$N $IT mesh/nghost=3 mesh/nx1=$N meshblock/nx1=$N \
        rad_m1/kappa_s=$K rad_m1/reconstruct=$rec time/tlim=0.5 output1/dt=0.05 &
  done
done
wait
echo RECON_DONE
