#!/bin/bash -l
# MILESTONE 3a2 gates.  Every command that produced a number in RESULTS.txt.  Each case
# is run, analysed and its dumps deleted immediately (viper inode quota).  TWO binaries,
# ONE AT A TIME (delete the other build directory first):
#   build_cpu_m1  : cmake -S . -B build_cpu_m1  -D CMAKE_BUILD_TYPE=Release
#   build_cpu_box : cmake -S . -B build_cpu_box -D CMAKE_BUILD_TYPE=Release \
#                         -D PROBLEM=box_convection
# with `module purge; module load gcc/14 cmake/4.0`, serial CPU.
set -u
R=/viper/u2/jinma/ATHENAK/bench/wt_m1impl2
T=$R/tests_m1
X=$R/build_cpu_m1/src/athena
XB=$R/build_cpu_box/src/athena
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IJ=$R/inputs/tests/rad_m1_jump.athinput
IM=$R/inputs/tests/rad_m1_marshak.athinput
I1D=$R/inputs/tests/rad_m1_pulse1d.athinput
IA=$R/inputs/tests/rad_m1_atmosphere.athinput
cd $T

# ---------------------------------------------------------------- STEP 0, explicit path
#   $X -i $IT rad_m1/kappa_s=128000 rad_m1/reconstruct=plm time/tlim=12.0 output1/dt=1.2
#   python3 t3_pulse.py <d>/bin/*.bin --c 1 --rho 1 --kappa 128000 --subtract-min --quiet
#       -> 1.000015, identical to runs_1cB and runs_3a.

# ------------------------------------------------------------------- I1, the thick pulse
#   $X -i $IT rad_m1/kappa_s={1280,128000,128000000} rad_m1/transport=implicit_x1 \
#       rad_m1/implicit_cfl=C time/cfl_number=1e6 rad_m1/implicit_flux={central,berthon,ap_hll}
#   python3 t3_pulse.py <d>/bin/*.bin --c 1 --rho 1 --kappa <k> --subtract-min --quiet

# --------------------------------------------------------- I6, the free-streaming pulse
#   $X -i $I1D rad_m1/transport=implicit_x1 rad_m1/implicit_cfl=C time/cfl_number=1e6 \
#       rad_m1/implicit_flux={central,berthon} [ rad_m1/implicit_recon=plm_dc ]
# the amplitude is (max-min)(final)/(max-min)(initial) of column 4 of the two tab files.

# --------------------------------------------------------------- I2, the opacity jump
T3B="--c 1 --kappa 0.64 --ratio 1000 --flux 1.0e-3 --e-left 1.0 --quiet"
IMPJ="rad_m1/transport=implicit_x1 rad_m1/implicit_bc_x1min=efix \
      rad_m1/implicit_bc_x1max=efix time/cfl_number=1e6"
#   $X -i $IJ $IMPJ rad_m1/implicit_cfl={1,100} rad_m1/implicit_flux=F
#   python3 t3b_jump.py <d>/tab/*.m1.00004.tab $T3B

# ------------------------------------------------------------------ I4, the Marshak wave
T6="--table $T/runs_3a/t6_ref_sn.txt --marshak --a-rad 1e30 --cv 1.5 --rho 1 \
    --centre 0.0 --quiet"
IMPM="rad_m1/transport=implicit_x1 rad_m1/implicit_bc_x1min=marshak \
      rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak \
      time/cfl_number=1e6"
#   $X -i $IM [ $IMPM rad_m1/implicit_cfl=C rad_m1/implicit_flux=F rad_m1/implicit_recon=RC ]
#   python3 t6_marshak.py <d>/tab/*.m1.00001.tab --hydro <d>/tab/*.hydro_w.00001.tab $T6

# ------------------------------------------------- I6b, the static grey atmosphere (NEW)
#   $X -i $IA rad_m1/transport=implicit_x1 rad_m1/implicit_cfl={1,100,1e4} \
#       time/cfl_number=1e6 rad_m1/implicit_flux=F time/tlim=2000 output1/dt=2000 \
#       output2/dt=2000
#   python3 t9_atmosphere.py <d>/tab/*.m1.00001.tab [--json plots/impl2_atm_F.json]
#   EXPLICIT reference on a thinner column (tau_bot ~ 10; the explicit scheme needs
#   ~tau^2 light crossings):
#   $X -i $IA problem/atm_scale_h=0.162 time/tlim=300 output1/dt=300 output2/dt=300
#   python3 t9_atmosphere.py <d>/tab/*.m1.00001.tab --scale-h 0.162
#   steps to steady state: re-run with time/nlim=N, N = 1,2,4,... and compare the last
#   tab of consecutive N until max|dE/E| < 1e-10.

# --------------------------------------------------------- I8 / LIMIT 3, the He column
# needs build_cpu_box:
#   $XB -i $T/runs_3a/he_box_m1_1d_impl.athinput problem/m1_top_bc=dark \
#       output2/dt=1.0e30 output3/dt=1.0e30 time/tlim=1000.0 \
#       rad_m1/implicit_bmom_half={false,true}
#   python3 runs_3a2/he_v1.py <d>/m1col.user.hst [--json plots/impl2_he_*.json]
echo "see RESULTS.txt; this file documents the commands, it does not run them"
