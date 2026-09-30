# sourced: arm definitions for the m1-positivity decisive runs (user brief 09-30)
W=/viper/ptmp2/jinma/m1pos_0930
BIN=$W/bin/athena_gpu_he_star_m1_7192a176     # m1-positivity 7192a176
K128="mesh/nx2=128 mesh/nx3=128 meshblock/nx2=64 meshblock/nx3=64 rad_m1/implicit_opac_newton=true problem/mlt_closure=adaptive problem/mlt_relax_time=1500 problem/mlt_closure_dout=587.5 problem/mlt_ramp_start=23500 problem/mlt_ramp_time=4700"
K64="rad_m1/implicit_opac_newton=true problem/mlt_closure=adaptive problem/mlt_relax_time=1500 problem/mlt_closure_dout=587.5 problem/mlt_ramp_start=23500 problem/mlt_ramp_time=4700"
OUT="output2/dt=2350 output3/dt=1175 output4/dt=1175"
G6="rad_m1/implicit_opac_newton_guard_mode=6"
declare -A R KF KX
# 128^2 (ad3d_128 rst 11, t = 25850)
R[off128]=$W/rst11.rst;  KF[off128]=off; KX[off128]="$K128"
R[onv128]=$W/rst11.rst;  KF[onv128]=on;  KX[onv128]="$K128 $G6"
R[onn128]=$W/rst11.rst;  KF[onn128]=on;  KX[onn128]="$K128 $G6 rad_m1/implicit_vimp=false"
R[be128]=$W/rst11.rst;   KF[be128]=off;  KX[be128]="$K128 rad_m1/time_scheme=be rad_m1/implicit_vimp=false"
# 64^2 hllc (ad3d_hllc rst 12, t = 28200)
R[offh64]=$W/rsth12.rst; KF[offh64]=off; KX[offh64]="$K64"
R[onnh64]=$W/rsth12.rst; KF[onnh64]=on;  KX[onnh64]="$K64 $G6 rad_m1/implicit_vimp=false"
R[onvh64]=$W/rsth12.rst; KF[onvh64]=on;  KX[onvh64]="$K64 $G6"
# run_arm <arm> <dir> <tlim> <walllimit> [extra keys]
run_arm() { local a=$1 d=$2 tl=$3 wl=$4; shift 4
  mkdir -p $d; cd $d
  local r=${R[$a]}; local last=$(ls -t $d/rst/*.rst 2>/dev/null | head -1); [ -n "$last" ] && r=$last
  echo "#### $a from $r tlim=$tl $(date +%T)"
  srun -n 2 $BIN -r $r -i $W/keys/${KF[$a]}.athinput -d $d -t $wl ${KX[$a]} time/tlim=$tl "$@" >> $d/run.log 2>&1
  echo "#### $a rc=$? $(date +%T)"; }
R[beon128]=$W/rst11.rst; KF[beon128]=on; KX[beon128]="$K128 $G6 rad_m1/time_scheme=be rad_m1/implicit_vimp=false"
