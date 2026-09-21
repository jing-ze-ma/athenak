#!/bin/bash
# Restart each relaxation run from its LAST restart file, advance one cycle, and dump one
# radial column, for three columns: night (mu0 = -0.922), day (mu0 = +0.922) and day
# (mu0 = +0.382).  The correlated-k column dump is one-shot and one column per run, which
# is why this is six short restarts rather than one.
set -e
cd "$(dirname "$0")"
BIN=/viper/u2/jinma/ATHENAK/athenak/build_cksph_new/src/athena
for sw in off on; do
  RUN=g_relax_$sw
  # the clean-exit restart, i.e. the final state.  The <output4> rst cadence in the input
  # is 3.05e6 s, far beyond this run, so the only restarts are the t = 0 one and the one
  # written at exit; the two runs stop at slightly different simulated times (both ran
  # 4000 cycles, and their dt histories differ by a few per cent), which is reported.
  RST=$(ls -1 $RUN/rst/*.rst | tail -1)
  for col in "night 0 2" "day92 1 2" "day38 0 4"; do
    set -- $col
    D=h_${1}_$sw
    rm -rf $D; mkdir -p $D
    # nlim is compared against the ABSOLUTE cycle count and the restart carries
    # nlim = 4000 = the cycle it stopped at, so it has to be raised; -t then stops the run
    # after a few seconds of wall clock.  The column dump is one-shot and fires on the
    # FIRST RT call, so it records the restarted state.
    ( cd $D && $BIN -r ../$RST -t 00:00:06 time/nlim=999999 \
        problem/ck_dump_file=col.txt problem/ck_dump_m=$2 problem/ck_dump_k=$3 \
        output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 > run.log 2>&1 ) \
      || { echo "FAIL $D"; tail -15 $D/run.log; exit 1; }
    rm -rf $D/bin $D/rst
    echo "$D  $(head -2 $D/col.txt | tail -1)"
  done
done
