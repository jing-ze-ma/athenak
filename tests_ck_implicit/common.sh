REF=/viper/u2/jinma/ATHENAK/bench/wt_ckimpl_ref/build_cpu/src/athena
NEW=/viper/u2/jinma/ATHENAK/bench/wt_ckimpl/build_cpu/src/athena
IN=/viper/u2/jinma/ATHENAK/bench/wt_ckimpl/inputs/tests/dhj_ck_implicit.athinput
BASE="output3/dt=1e30 output4/dt=1e30 problem/rt_use_cons=true"
run() { local d=$1; shift; local b=$1; shift
  rm -rf "$d"; mkdir -p "$d"; ( cd "$d" && $b -i $IN $BASE "$@" > run.log 2>&1 ); echo "$? $d"; }
