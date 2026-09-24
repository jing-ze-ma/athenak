#!/bin/bash
cd /viper/ptmp2/jinma/thinstab_0924
./run.sh off_s1 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni.athinput &
for ks in 0.25 0.5 1 2 4; do
  ./run.sh ctr15_k$ks s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/kappa_s=$ks &
done
./run.sh ctr15_k0.05 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/kappa_s=0.05 &
./run.sh ctr15_nyq s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 problem/atm_seed_k=0 &
./run.sh ctr15_ker s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/closure=kershaw &
./run.sh ctr15_q09_k2 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/marshak_q=0.9 rad_m1/kappa_s=2 &
./run.sh ctr15_q09_k025 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/marshak_q=0.9 rad_m1/kappa_s=0.25 &
./run.sh ctr10_k1 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.0 &
./run.sh ctr15_nu30 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/implicit_cfl=30 &
./run.sh ctr15_odn s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/implicit_offdiag=none &
./run.sh ctr15_k1_s1em2 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 problem/atm_seed=1.0e-2 &
./run.sh ref500_k1 s1 /viper/ptmp2/jinma/thinstab_0924/inp/uni_ctr.athinput time/nlim=500 output3/dcycle=25 &
wait
