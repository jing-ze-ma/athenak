#!/bin/bash
# Gate (c), the star-off L(r) flat test: restart each relaxed column from its last
# restart file and dump the NIGHT column (mu0 = -0.922, m = 0, k = 2) once, then
#   python3 ../tests_ck_sph/lprof.py h_night_f0/col.txt h_night_f1/col.txt h_night_f2/col.txt
# In radiative equilibrium L(r) = A(r) F_net(r) must be flat; the plane-parallel form
# fails it by exactly the area ratio (tests_ck_sph/README.md gate 3).
set -e
T=$(cd "$(dirname "$0")" && pwd)
BIN=$T/../build_cksd/src/athena
for f in 0 1 2; do
  RST=$(ls -1 $T/g_relax_f$f/rst/*.rst | tail -1)
  d=$T/h_night_f$f
  rm -rf $d && mkdir -p $d
  ( cd $d && $BIN -r $RST -t 00:00:20 time/nlim=999999999 \
      problem/ck_sweep_form=$f \
      problem/ck_dump_file=col.txt problem/ck_dump_m=0 problem/ck_dump_k=2 \
      output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 \
      > run.log 2>&1 ) || echo "FAIL h_night_f$f"
done
echo dumps done
