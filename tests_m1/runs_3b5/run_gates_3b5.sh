#!/bin/bash -l
# Milestone 3b phase D gates (docs/dev/rad_m1_implicit_design.md sect. 13).  Every case
# below was run exactly as written, from tests_m1/runs_3b5, on serial CPU; the dumps of a
# case are deleted as soon as it is measured (viper inode quota).
module purge; module load gcc/14 cmake/4.0
A=/viper/u2/jinma/ATHENAK/athenak
M1=$A/build_cpu_m1/src/athena          # default pgen: the thick-pulse gates
BOX=$A/build_cpu_box/src/athena        # -D PROBLEM=box_convection: the He slab
H=$A/tests_m1/runs_3b5
P="mesh/nx1=64 meshblock/nx1=64 mesh/nx2=64 meshblock/nx2=64 mesh/x2min=0.0 \
mesh/x2max=1.0 problem/pulse_y0=0.5 time/cfl_number=1e6 rad_m1/transport=implicit \
rad_m1/implicit_solver=bicgstab"

# ---- G-oblique: tau_cell = 0.1, the flux direction sweeps every angle -----------------
cd $H/ob_explicit && $M1 -i $H/pulse_md.athinput -d . $P rad_m1/transport=explicit \
  rad_m1/kappa_s=6.4 time/tlim=0.2 output1/dt=0.2 > run.log 2>&1
for m in none lagged operator; do
  cd $H/ob_$m && $M1 -i $H/pulse_md.athinput -d . $P rad_m1/kappa_s=6.4 \
    rad_m1/implicit_cfl=1.0 rad_m1/implicit_offdiag=$m time/tlim=0.2 output1/dt=0.2 \
    > run.log 2>&1
done
cd $H/ob_operator_step && $M1 -i $H/pulse_md.athinput -d . $P rad_m1/kappa_s=6.4 \
  rad_m1/implicit_cfl=1.0 rad_m1/implicit_offdiag=operator \
  rad_m1/implicit_closure_lag=step time/tlim=0.2 output1/dt=0.2 > run.log 2>&1
for d in ob_explicit ob_none ob_lagged ob_operator ob_operator_step; do
  python3 $H/oblique.py --label $d --dir $H/$d --ref $H/ob_explicit
done

# ---- G-pulse: the phase-B/C diffusion rate and isotropy -------------------------------
for m in lagged operator; do
  cd $H/p_tau10_$m && $M1 -i $H/pulse_md.athinput -d . $P rad_m1/kappa_s=640 \
    rad_m1/implicit_cfl=1.0 time/tlim=0.1984 output1/dt=0.0165333 \
    rad_m1/implicit_offdiag=$m > run.log 2>&1
  cd $H/p_tau1e3_$m && $M1 -i $H/pulse_md.athinput -d . $P rad_m1/kappa_s=64000 \
    rad_m1/implicit_cfl=1e4 time/tlim=1875.0 output1/dt=156.25 \
    rad_m1/implicit_offdiag=$m > run.log 2>&1
done
# ...and the same two with rad_m1/implicit_closure_lag=step (dirs p_*_opstep), plus the
# 3-D 19-point stencil (dirs p3d_lagged, p3d_operator, input pulse_md3.athinput):
#   mesh/nx{1,2,3}=32 rad_m1/kappa_s=32000 rad_m1/implicit_cfl=1.0 time/tlim=9.92
python3 $H/../runs_3b3/pulse_md.py $H/p_tau10_operator/bin/*.bin --kappa 640 --label x

# ---- G-static and G-seed: the 2-D He slab --------------------------------------------
S="$H/he_slab_m1_2d.athinput"
cd $H/stat_op_pass && $BOX -i $S -d . time/tlim=100.0 time/ndiag=200 \
  rad_m1/implicit_offdiag=operator > log.txt 2>&1
cd $H/stat_op_step && $BOX -i $S -d . time/tlim=100.0 time/ndiag=200 \
  rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step > log.txt 2>&1
cd $H/stat_none_step && $BOX -i $S -d . time/tlim=100.0 time/ndiag=200 \
  rad_m1/implicit_offdiag=none rad_m1/implicit_closure_lag=step > log.txt 2>&1
for m in operator none lagged; do
  cd $H/seed_$m && $BOX -i $S -d . time/tlim=200.0 problem/vpert=1.0e-3 time/ndiag=200 \
    rad_m1/implicit_offdiag=$m rad_m1/implicit_closure_lag=step > log.txt 2>&1
done
for d in stat_op_pass stat_op_step stat_none_step seed_operator seed_none; do
  cp $H/$d/log.txt $H/$d/run.log
done
python3 $H/../runs_3b4/heslab.py seeded $H/seed_op static $H/stat_op_step \
  --out $H/../plots/impl3b5_heslab.json

# ---- regression and style ------------------------------------------------------------
bash $A/tests_gate_merge/postmerge.sh
cd $A/tst/test_suite/style && ./check_athena_cpp_style.sh
