#!/bin/bash
# collect_arm.sh ARM [OHMIC_DT]
# Prints the stability / dt / wall / ms-per-cycle lines of one w121_mhd_0929 arm for RESULTS_deltaai_w121mhd.
# Reads runs/ARM/run.*.log (all restart links, in job-id order) and runs/ARM/ana_eta.txt if present.
# Read-only. OHMIC_DT (s, optional, e.g. 1.755 for 1e13 or 0.18 for 1e14) -> fraction of dt samples within 3 %
# of it (= Ohmic-limited); without it, "flat" = fraction of consecutive 100-cycle samples with identical dt
# (Ohmic dt is piecewise constant, the CFL dt is not).
set -u
R=/work/nvme/bivj/jma20/w121_mhd_0929
ARM=${1:?usage: collect_arm.sh ARM [OHMIC_DT]}
OHM=${2:-0}
D=$R/runs/$ARM
[ -d "$D" ] || { echo "no arm dir $D"; exit 1; }
T0=3.304605e7; PROT=110153.5
# tlim is set on the athena command line by run_arm_deltaai.sub (the input's tlim is rot 300): take it from the job .out
TLIM=$(grep -ho 'time/tlim=[0-9.e+]*' $R/logs/w${ARM}.*.out 2>/dev/null | tail -1 | cut -d= -f2)
[ -n "$TLIM" ] || TLIM=33266357.0   # rot 302 = T0 + 2 P_rot
LOGS=$(ls $D/run.*.log 2>/dev/null | sort -t. -k2,2n)
[ -n "$LOGS" ] || { echo "$ARM: no run.*.log (not run on DeltaAI)"; exit 0; }
echo "== $ARM  logs: $(for f in $LOGS; do basename $f; done | tr '\n' ' ')"
for f in $LOGS; do
  echo "   $(basename $f): ranks $(grep -m1 -o 'Number of parallel ranks = [0-9]*' $f | awk '{print $NF}'), end: $(grep -m1 -E 'Terminating on|FATAL' $f || echo 'RUNNING/killed (no Terminating line)')"
done
cat $LOGS | awk '
  /NOT-CONVERGED/ {nc++} /FATAL/ {fa++} tolower($0) ~ /(^|[^a-z])nan([^a-z]|$)/ {nan++} tolower($0) ~ /dt collapse/ {dc++}
  END {printf "   stability: NOT-CONVERGED %d  FATAL %d  nan %d  dtCOLLAPSE %d\n", nc, fa, nan, dc}'
OUT=$(for f in $LOGS; do
  awk -v f=$(basename $f) -v ohm=$OHM -v P=$PROT '
    /^elapsed=[^ ]+ cycle=[0-9]+ time=[^ ]+ dt=/ {
      split($0, a, /[= ]/); el=a[2]+0; cy=a[4]+0; t=a[6]+0; dt=a[8]+0
      if (n==0) {el0=el; cy0=cy; t0=t; dmin=dt; dmax=dt; dfirst=dt}
      n++; s+=dt; if (dt<dmin) dmin=dt; if (dt>dmax) dmax=dt
      if (n>1 && dt==dprev) flat++
      if (ohm>0 && dt>0.97*ohm && dt<1.03*ohm) nohm++
      dprev=dt; el1=el; cy1=cy; t1=t; E[n]=el; C[n]=cy; T[n]=t; DT[n]=dt
    }
    /RKG super-stepping total stage=/ {split($0,b,/[= ]/); st[b[5]+0]++; ns++; dd+=b[7]}
    /cpu time used/ {cpu=$NF}
    END {
      if (n<2) {print "   " f ": <2 cycle lines"; exit}
      printf "   %s: cycles %d-%d  t %.6e-%.6e  dt first %.4f last %.4f mean %.4f min %.4f max %.4f  flat %.2f", f, cy0, cy1, t0, t1, dfirst, dprev, s/n, dmin, dmax, (flat+0)/(n-1)
      if (ohm>0) printf "  f(dt~Ohmic %.3g) %.2f", ohm, (nohm+0)/n
      printf "\n   %s: wall(last elapsed) %.1f s  cpu %s s  ms/cycle %.2f (excl. startup)\n", f, el1, (cpu==""?"-":cpu), 1e3*(el1-el0)/(cy1-cy0)
      for (i=1; i<=n && T[i] < t1-0.5*P; i++); if (i<n) { sd=0; for (j=i; j<=n; j++) sd+=DT[j];
        printf "   %s: last 0.5 rot (t>=%.4e): %.1f min/rot, %.2f ms/cycle, dt mean %.4f s\n", f, T[i], (E[n]-E[i])/60/((T[n]-T[i])/P), 1e3*(E[n]-E[i])/(C[n]-C[i]), sd/(n-i+1) }
      if (ns>0) {printf "   %s: RKG stages s:", f; for (k in st) printf " %d:%.3f", k, st[k]/ns; printf "  (n=%d, mean dt_diff %.4f)\n", ns, dd/ns}
      printf "SEG %.6f %.6f %d %d %.6f\n", el1, (t1-t0), cy1-cy0, cy1, t1
    }' $f
done)
echo "$OUT" | grep -v '^SEG'
echo "$OUT" | awk -v T0=$T0 -v P=$PROT -v TL=$TLIM '
  /^SEG/ {w+=$2; dt+=$3; dc+=$4; tend=$6}
  END {if (dt>0) printf "   TOTAL: %.3f rot of %.3f (%.1f %%), wall %.1f s -> %.1f min/rot, %.2f ms/cycle (wall incl. startup / cycles), mean dt (time/cycles) %.4f s\n",
        (tend-T0)/P, (TL-T0)/P, 100*(tend-T0)/(TL-T0), w, w/60/(dt/P), 1e3*w/dc, dt/dc}'
if [ -f $D/ana_eta.txt ]; then
  echo "   ana_eta.txt summary:"; grep -E '^# (arm|cost)|^## |jet max|E_mag|div B' $D/ana_eta.txt | sed 's/^/     /'
else
  echo "   ana_eta.txt: none (run: PYTHONPATH=<dir with empty h5py.py> python3 deltaai_pkg/ana_eta.py runs/$ARM)"
fi
