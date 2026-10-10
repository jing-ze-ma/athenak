#!/bin/bash
# xthinfix-1009 STEEP battery (Caltech, 10-10) (correctness only; NO timings from here).
# usage: REPO=<athenak checkout at xthinfix-1009> ATHENA_CPU=<none CPU binary> XTF_RUN=<run dir> P=<cores> \
#        bash battery_cpu.sh [order|gates|beams|ana|all]
# Arms: cen central; hr production (xthin 30, mode all); hrx0 xthin 0 (reference);
#       hrs / hrs15 / hrs60 steep gate X0 30 / 15 / 60 (dc); hrsp steep + plm; hrp all + plm.
set -u
B=$(cd "$(dirname "$0")" && pwd)
: "${REPO:?}" "${ATHENA_CPU:?}" "${XTF_RUN:?}"
P=${P:-16}
export REPO ATHENA_CPU XTF_RUN
mkdir -p "$XTF_RUN"
WHAT=${1:-all}
CEN="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
HRK="rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange"
IDF="rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6"
declare -A A
A[cen]="rad_m1/implicit_flux=central"
A[hr]="$HRK $IDF rad_m1/implicit_blend_xthin=30"
A[hrb]="$HRK $IDF rad_m1/implicit_blend_xthin=30 rad_m1/implicit_blend_xthin_mode=beam"
A[hrx0]="$HRK $IDF rad_m1/implicit_blend_xthin=0"
ST="$HRK $IDF rad_m1/implicit_blend_xthin_mode=steep"
A[hrs]="$ST rad_m1/implicit_blend_xthin=30"
A[hrs15]="$ST rad_m1/implicit_blend_xthin=15"
A[hrs60]="$ST rad_m1/implicit_blend_xthin=60"
A[hrsp]="$ST rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm"
A[hrp]="${A[hr]} rad_m1/implicit_hr_recon=plm"

# ---- 1. convergence ladders (od.py): pulses tau_cell 2.5e-4..400, grey atmosphere vs Hopf 32-512,
#         coupled radiation-acoustic wave P 100 tau_lambda 10 / 1e3; be and hesdirk2; X (space), T (time), C (combined)
order() {
  local J=$XTF_RUN/jobs_order.txt; : > $J
  python3 $B/ana/od.py jobs atm:cen,hr,hrx0,hrs,hrs15,hrs60,hrsp,hrp:be:S:1 >> $J
  for c in pulse_k0.128 pulse_k12.8 pulse_k128 pulse_k1280 pulse_k12800 rw_t10 rw_t1000; do
    python3 $B/ana/od.py jobs $c:cen,hr,hrs,hrs15,hrs60,hrsp:be,hesdirk2:XTC:1 >> $J
  done
  python3 $B/ana/od.py run $J $P > $XTF_RUN/logs_order.txt 2>&1
}

# ---- 2. G1 (grey atmosphere vs Hopf, cfl 1e2/1e4), G5 (thick pulse diffusion rate), G3 (Marshak)
gates() {
  local G=$XTF_RUN/gates; mkdir -p $G
  python3 - "$REPO" "$G" <<'PY'
import sys
repo, g = sys.argv[1], sys.argv[2]
keys = ('implicit_flux_faces   = x1\nimplicit_flux_beam    = closure\nimplicit_blend_alpha  = 1.0\n'
        'implicit_blend_r0     = 1.5\nimplicit_blend_xthin = 0.0\nimplicit_blend_xthin_mode = all\n'
        'implicit_blend_xthin_wmin = 0.0\nimplicit_hr_recon = dc\nimplicit_hr_recon_fresh = 2\nimplicit_hr_damp = 0.0\nvet_nmu      = 8\nvet_nphi     = 32\nvet_source_noesrc = false')
s = open(repo + '/inputs/tests/rad_m1_thick_pulse.athinput').read()
s = s.replace('nx2       = 1\nx2min     = -0.5', 'nx2       = 4\nx2min     = -0.5', 1)
a = s.index('<meshblock>')
s = s[:a] + s[a:].replace('nx2       = 1', 'nx2       = 4', 1)
s = s.replace('implicit_flux         = central', 'implicit_flux         = central\n' + keys, 1)
s = s.replace('implicit_solver       = line_jacobi', 'implicit_solver       = bicgstab', 1)
s = s[:s.index('<output2>')]
open(g + '/pulse2d.athinput', 'w').write(s)
s = open(repo + '/inputs/tests/rad_m1_marshak.athinput').read()
s = s.replace('nx2       = 1', 'nx2       = 4', 2)
s = s.replace('implicit_flux         = central', 'implicit_flux         = central\n' + keys, 1)
s = s.replace('variable    = m1\n', 'variable    = m1\nslice_x2    = 0.0\n', 1)
s = s.replace('variable    = hydro_w\n', 'variable    = hydro_w\nslice_x2    = 0.0\n', 1)
open(g + '/marsh2d.athinput', 'w').write(s)
PY
  local T9=$REPO/tests_m1/t9_atmosphere.py T3=$REPO/tests_m1/t3_pulse.py T6=$REPO/tests_m1/t6_marshak.py
  local REF=$REPO/tests_m1/runs_3a/t6_ref_sn.txt
  local L=$G/jobs.txt; : > $L
  for a in hr hrs hrs15 hrs60 hrsp; do
    for c in 1e2 1e4; do echo "g1|$G/g1_${a}_$c|$B/inp/atm2d.athinput|$CEN ${A[$a]} rad_m1/implicit_cfl=$c time/tlim=2000 output1/dt=2000 output2/dt=2000" >> $L; done
    for k in 12.8 128 1280 12800; do echo "g5|$G/g5_${a}_k$k|$G/pulse2d.athinput|$CEN ${A[$a]} rad_m1/implicit_cfl=1 rad_m1/kappa_s=$k" >> $L; done
    for c in 1 10; do echo "g3|$G/g3_${a}_c$c|$G/marsh2d.athinput|$CEN ${A[$a]} rad_m1/implicit_cfl=$c rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak" >> $L; done
  done
  export T9 T3 T6 REF B
  xargs -P $P -d '\n' -n 1 bash -c '
    IFS="|" read -r g R IN K <<< "$0"; mkdir -p $R; cd $R || exit 1
    OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
    rc=$?; pic=$(grep -o "Picard iterations mean=[0-9.e+-]*" run.log | tail -1)
    case $g in
      g1) r=$(python3 $B/ana/t9h.py $(ls tab/*.m1.*.tab | tail -1) 2>&1) ;;
      g5) k=${R##*_k}; r=$(python3 $T3 bin/*.bin --c 1 --rho 1 --kappa $k --subtract-min --quiet --reduce mean 2>&1) ;;
      g3) r=$(python3 $T6 $(ls tab/*.m1.*.tab | tail -1) --hydro $(ls tab/*.hydro_w.*.tab | tail -1) --table $REF --marshak --a-rad 1e30 --cv 1.5 --rho 1 --centre 0.0 --quiet 2>&1) ;;
    esac
    echo "$(basename $R) rc=$rc $pic | $(echo $r | tr "\n" " ")"' < $L > $G/gates.out 2>&1
}

# ---- 3. beam stability / accuracy: xb20 ba0 ba20 cyl shd3b, c dt/dx ~300, 100 cycles, np 1; arms hr hrs hrs15 hrs60 hrsp hrp (+cen)
beams() {
  local D=$XTF_RUN/beams L=$XTF_RUN/beams/jobs.txt; mkdir -p $D; : > $L
  for t in shd3b cyl xb20 ba0 ba20; do
    for a in hr hrs hrs15 hrs60 hrsp hrp cen; do
      K="$CEN"; [ $a != cen ] && K="$CEN ${A[$a]}"
      echo "$D/$t/$a|$B/inp/$t.athinput|time/nlim=100 rad_m1/c_light=1000 $K" >> $L
    done
  done
  xargs -P $P -d '\n' -n 1 bash -c '
    IFS="|" read -r R IN K <<< "$0"; mkdir -p $R; cd $R || exit 1
    /usr/bin/time -f "%e s wall" -o time.txt env OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R $K > run.log 2>&1
    echo "rc=$? fatal=$(grep -c FATAL run.log)" > rc.txt' < $L
  cp $B/ana/exact_shd3b.npz $D/shd3b/
}

ana() {
  local O=$XTF_RUN/RESULTS; mkdir -p $O
  { for st in atm:cen,hr,hrx0,hrs,hrs15,hrs60,hrsp,hrp:be:S:1; do python3 $B/ana/od.py eval atm:cen,hr,hrx0,hrs,hrs15,hrs60,hrsp,hrp:be:S:1; done
    for c in pulse_k0.128 pulse_k12.8 pulse_k128 pulse_k1280 pulse_k12800; do
      python3 $B/ana/od.py eval $c:cen,hr,hrs,hrs15,hrs60,hrsp:be,hesdirk2:XTC:1; done
    python3 $B/ana/od.py eval rw_t10:cen,hr,hrs,hrs15,hrs60,hrsp:be,hesdirk2:XTC:1
    python3 $B/ana/od.py eval rw_t1000:cen,hr,hrs,hrs15,hrs60,hrsp:be,hesdirk2:XTC:1 "X=32,64,128,256;C=32,64,128,256;T=2,4,8,16"
  } > $O/order_eval.txt 2>&1
  python3 $B/ana/ana_atm.py cen hr hrx0 hrs hrs15 hrs60 hrsp hrp > $O/atm_hopf.txt 2>&1
  python3 $B/ana/ana_hrcen.py hr hrs hrs15 hrs60 hrsp > $O/arm_minus_cen.txt 2>&1
  cp $XTF_RUN/gates/gates.out $O/gates.txt 2>/dev/null
  local D=$XTF_RUN/beams
  { for t in shd3b cyl xb20 ba0 ba20; do
      for a in hr hrs hrs15 hrs60 hrsp hrp cen; do R=$D/$t/$a
        echo "$t $a $(cat $R/rc.txt) $(cat $R/time.txt) | $(grep -o 'Picard iterations mean=[0-9.e+-]* max=[0-9.e+-]* NON-CONVERGED=[0-9.e+-]*' $R/run.log | tail -1) | $(grep -o 'inner iterations mean=[0-9.e+-]* max=[0-9.e+-]*' $R/run.log | tail -1) | $(grep -o 'breakdowns=[0-9.e+-]*' $R/run.log | tail -1) | $(grep -o 'plm positivity fallbacks=[0-9.e+-]*' $R/run.log | tail -1)"
      done; done
    python3 $B/ana/cylref.py $B/inp/cyl.athinput $D/cyl/hr $D/cyl/hrs $D/cyl/hrsp $D/cyl/cen 2>&1 | grep -v '^    \|y: '
    for t in xb20 ba0 ba20; do python3 $B/ana/collref.py $B/inp/$t.athinput $D/$t/hr $D/$t/hrs $D/$t/hrsp $D/$t/cen 2>&1; done
    python3 $B/ana/shdfs.py $D/shd3b/exact_shd3b.npz $D/shd3b/hr $D/shd3b/hrs $D/shd3b/hrsp $D/shd3b/cen 2>&1
  } > $O/beams.txt
  echo "results in $O"
}

case $WHAT in
  order) order ;; gates) gates ;; beams) beams ;; ana) ana ;;
  all) order; gates; beams; ana ;;
esac
