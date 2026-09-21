#!/bin/bash -l
# MILESTONE 3a gates: <rad_m1>/transport = implicit_x1.  Every command that produced a
# number in RESULTS.txt is here.  Each case is run, analysed and its dumps deleted
# immediately: the viper inode quota is at its limit, and a parallel batch of these runs
# hits "Disk quota exceeded" in the middle of an output.
#
# Two binaries are needed, ONE AT A TIME (delete the build directory before the other):
#   build_cpu_m1  : cmake -S . -B build_cpu_m1  -D CMAKE_BUILD_TYPE=Release
#   build_cpu_box : cmake -S . -B build_cpu_box -D CMAKE_BUILD_TYPE=Release \
#                         -D PROBLEM=box_convection
# with `module purge; module load gcc/14 cmake/4.0`, serial CPU.
#
# NOTE on the timestep.  In implicit mode the module never sub-cycles: it takes the mesh
# step.  <rad_m1>/implicit_cfl > 0 makes the module own the mesh step and asks for a
# RADIATION CFL c dt/dx; because Mesh::NewTimeStep takes the MINIMUM over the modules,
# the tests that want a large radiation CFL also pass time/cfl_number=1e6 so that the
# hydro step (a static, uniform medium in every one of them) cannot limit it.
set -u
R=/viper/u2/jinma/ATHENAK/bench/wt_m1impl
T=$R/tests_m1
X=$R/build_cpu_m1/src/athena
XB=$R/build_cpu_box/src/athena
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IJ=$R/inputs/tests/rad_m1_jump.athinput
IM=$R/inputs/tests/rad_m1_marshak.athinput
IE=$R/inputs/tests/rad_m1_equil.athinput
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
IU=$R/inputs/tests/rad_m1_advect_uniform.athinput
I1D=$R/inputs/tests/rad_m1_pulse1d.athinput
cd $T

# ---------------------------------------------------------------- step 0: the EXPLICIT
# scheme must be unchanged.  runs_3a/step0.sh, then:
#   python3 t3_pulse.py runs_3a/step0/t3_tau1e3_plm/bin/*.bin --c 1 --rho 1 \
#       --kappa 128000 --subtract-min --quiet
#   A5="--c 2.99792458e10 --rho 1.2 --kappa 500 --v 3.0e7 --drift-tol 1 --quiet"
#   python3 t4_advect.py <t4_d512_split>/tab/*.m1.*.tab --hydro <...>/*.hydro_w.*.tab \
#       $A5 --static <t4_s512>/tab/m1_advect_pulse.m1.00004.tab
#   python3 t4_advect.py <t4b_full|t4b_ovc> ... --v 2.99792458e8 --kappa 83333.33333333333

# ------------------------------------------------------------------------ I1, run_i1.sh
# plus the four cases whose default tlim is shorter than one implicit step:
#   i1x_tau1e1_c100  kappa 1280   tlim 8     odt 0.8
#   i1x_tau1e3_c1e4  kappa 128000 tlim 1200  odt 120
#   i1x_nyq_c1e4     kappa 128000 tlim 780   odt 78    problem/nyquist_amp=1.0e-3
#   i1x_tau1e1_c1e4  kappa 1280   tlim 800   odt 80    (saturates, see RESULTS.txt)

# ---------------------------------------------------------------------------------- I2
# T3b opacity jump.  The implicit run holds the two END CELLS at their initial E
# (implicit_bc = efix), which is the face-flux form of the fixed-E ghost the explicit
# test uses and the only boundary pair that anchors the level of E in a pure-scattering
# column.
T3B="--c 1 --kappa 0.64 --ratio 1000 --flux 1.0e-3 --e-left 1.0 --quiet"
IMPJ="rad_m1/transport=implicit_x1 rad_m1/implicit_bc_x1min=efix \
      rad_m1/implicit_bc_x1max=efix time/cfl_number=1e6"
#   $X -i $IJ                                  # explicit, plm
#   $X -i $IJ rad_m1/reconstruct=dc            # explicit, dc
#   $X -i $IJ $IMPJ rad_m1/implicit_cfl={1,10,100}
#   python3 t3b_jump.py <dir>/tab/*.m1.00004.tab $T3B
# T3c top-hat: $X -i $IT problem/m1_test=tophat rad_m1/kappa_s=128000 time/tlim=12
#   output1/dt=1.2 [ rad_m1/transport=implicit_x1 rad_m1/implicit_cfl=C
#   time/cfl_number=1e6 ]; analysed with t3_pulse.py --tol 0.05

# ---------------------------------------------------------------------------------- I3
# equilibration, ideal gas (rad_m1_equil) and tabulated He EOS (rad_m1_equil_table).
# implicit_cfl plays exactly the role cfl_rad plays for the explicit scheme there:
# dt = implicit_cfl*dx/c = 33.356*implicit_cfl s, and the thermal coupling time is
# 1.4693e-8 s, so 4.4e-7 / 4.4e-5 / 4.4e-3 are 1e3 / 1e5 / 1e7 coupling times.
AN="--rho 1.0e-7 --kappa-p 0.4 --a-rad 2.0517256e-47 --erad0 1.0e12 --cv 1.5 --quiet"
#   $X -i $IE rad_m1/transport=implicit_x1 rad_m1/implicit_cfl=C problem/gas_eint={1e2,1e10}
#   python3 t5_equil.py <d>/tab/*.m1.*.tab --hydro <d>/tab/*.hydro_u.*.tab $AN --egas0 E0
# NOTE --cv 1.5: AthenaK's ideal gas has e = rho T/(gamma-1) in CODE temperature, so the
# script's c_v is 1/(gamma-1); with the default --mu 0.6 the reference is meaningless
# (it reports dE_eq = 9.1e15 for a run that IS at equilibrium -- the explicit scheme
# fails the same way, so this is a script invocation, not a scheme property).

# ---------------------------------------------------------------------------------- I4
T6="--table $T/runs_3a/t6_ref_sn.txt --marshak --a-rad 1e30 --cv 1.5 --rho 1 \
    --centre 0.0 --quiet"
IMPM="rad_m1/transport=implicit_x1 rad_m1/implicit_bc_x1min=marshak \
      rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak \
      time/cfl_number=1e6"
#   $X -i $IM mesh/nx1=N meshblock/nx1=N [ $IMPM rad_m1/implicit_cfl=C ]
#   python3 t6_marshak.py <d>/tab/*.m1.00001.tab --hydro <d>/tab/*.hydro_w.00001.tab $T6

# ------------------------------------------------------------------------------- I5, I7
TL="time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6"
N512="mesh/nx1=512 meshblock/nx1=512"
DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7"
#   $X -i $IP $TL $N512 $DYN rad_m1/transport=implicit_x1 [rad_m1/implicit_cfl=10
#       time/cfl_number=1e6]          # dt = dt_hydro when implicit_cfl is left at -1
#   python3 t4_advect.py ... $A5      # centre / width / drift
#   python3 t8_conserve.py <d>/tab/*.m1.*.tab --hydro <d>/tab/*.hydro_w.*.tab \
#       --c 2.99792458e10 --tol 1e-12 --quiet
#   $X -i $IU rad_m1/transport=implicit_x1                       # T4b
# restart (I7): a first run with time/nlim=95 and an <output3> rst block, then
#   $X -r <rst> -i <same input> <same overrides> time/tlim=4.8e-6
# (a command-line override of a block the RESTART file carries but the -i file does not
# is rejected, so the restart must be given the input file again).
# 1 vs 4 MeshBlocks along x2 (I7): mesh/nx2=16 with meshblock/nx2=16 vs 4 and
#   rad_m1/implicit_allow_multid=true; <mesh> needs >= 4 cells per block per dimension.
# free streaming (I6): $X -i $I1D mesh/nx1=128 meshblock/nx1=128 [implicit]; the damping
# is read off as the ratio of the final to the initial peak-to-trough amplitude.

# ---------------------------------------------------------------------------------- I8
# THE 1-D He COLUMN AT TRUE c.  Needs build_cpu_box.
#   $XB -i $T/runs_3a/he_box_m1_1d_impl.athinput problem/m1_top_bc=dark \
#       output2/dt=1.0e30 output3/dt=1.0e30 time/tlim=1000.0
#   control: the same input with rad_m1/transport=explicit rad_m1/chat_over_c=2.847e-3
#            rad_m1/subcycle=true      (= the milestone 2a G1 arm, K = 10)
# and the box G1 bitwise gate (build_cpu_box), tests_gate_merge/g1_gate.sh's run line:
#   $XB -i bench/hestar_fecz/box_w8/he_box_w8.athinput mesh/nx2=16 mesh/nx3=16 \
#       meshblock/nx1=134 meshblock/nx2=8 meshblock/nx3=8 time/nlim=50 \
#       [ problem/rt_implicit_column=0 problem/rt_col3_skip_sweep=false ]
#   cmp against athenak/tests_gate_merge/g1_new_m{3,0}
echo "see RESULTS.txt; this file documents the commands, it does not run them"
