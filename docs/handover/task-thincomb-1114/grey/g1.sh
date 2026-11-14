#!/bin/bash
# G1 static grey atmosphere (milestone 3c gate) with the hrup-1009 binary: central, berthon
# (closure HLL), half-range berthon, half-range blend tau (W1-hr), at implicit_cfl 1, 1e2,
# 1e4; tlim 2000; t9_atmosphere.py on the last tab dump.  Target <= 2e-3 (central 1.353e-3)
source /etc/profile.d/modules.sh; module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
D=/viper/ptmp2/jinma/blendall_1009/gates
X=/viper/ptmp2/jinma/blendall_1009/bc_hr_none/src/athena
T=/viper/ptmp2/jinma/wt_hrup_1009/tests_m1/t9_atmosphere.py
IN=$D/atm.athinput
python3 - <<'PY'
s = open('/viper/ptmp2/jinma/wt_hrup_1009/inputs/tests/rad_m1_atmosphere.athinput').read()
s = s.replace('implicit_flux         = central', 'implicit_flux         = central\nimplicit_flux_beam    = closure', 1)
open('/viper/ptmp2/jinma/blendall_1009/gates/atm.athinput', 'w').write(s)
PY
declare -A A
A[central]="rad_m1/implicit_flux=central"
A[berthon]="rad_m1/implicit_flux=berthon"
A[hr_be]="rad_m1/implicit_flux=berthon rad_m1/implicit_flux_beam=halfrange"
A[hr_w1]="rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau rad_m1/implicit_flux_beam=halfrange"
A[hr_be_step]="rad_m1/implicit_flux=berthon rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_recon_lag=step"
A[hr_w1_step]="rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_recon_lag=step"
for a in central berthon hr_be hr_w1 hr_be_step hr_w1_step; do
  for c in 1 1e2 1e4; do
    R=$D/g1/${a}_$c
    mkdir -p $R
    ( cd $R || exit 1
      nice -n 10 $X -i $IN -d $R rad_m1/transport=implicit_x1 time/cfl_number=1e6 ${A[$a]} \
        rad_m1/implicit_cfl=$c time/tlim=2000 output1/dt=2000 output2/dt=2000 > run.log 2>&1
      echo "$a cfl $c rc=$? $(grep -o 'Picard iterations mean=[0-9.e+-]*' run.log | tail -1) | $(python3 $T $(ls tab/*.m1.*.tab | tail -1) --quiet --label $a 2>&1 | tr '\n' ' ')" ) &
  done
  wait
done
