#!/bin/bash
cd /viper/ptmp2/jinma/s2_0924
M="rad_m1/closure=m1 rad_m1/implicit_precond=line"
./run.sh sym_m1d new none 4 $PWD/inp/sph_sym.athinput $M rad_m1/implicit_lin_tol=1.0e-14 rad_m1/implicit_tol=1.0e-12 &
./run.sh sym_m1a new none 4 $PWD/inp/sph_sym.athinput $M &
./run.sh sym_ker new none 4 $PWD/inp/sph_sym.athinput rad_m1/closure=kershaw rad_m1/implicit_precond=line &
wait
