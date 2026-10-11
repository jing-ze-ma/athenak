#!/bin/bash
# suite-1011 crooked pipe (JSD12 Sect 5.2, Cartesian), arms REC / REC-be, vet_sc.
# usage: BIN=<athena CPU binary> NT=<omp threads per run> bash run_pipe.sh <outdir>
# Runs the 4 cases below in parallel (one process each, OpenMP NT threads, 1 MeshBlock).
set -u
H=$(cd $(dirname $0) && pwd); OUT=${1:?outdir}; BIN=${BIN:?}; NT=${NT:-16}; PY=${PY:-python3}
CEN="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
HRK="rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange"
IDF="rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6"
TJE="$CEN $HRK $IDF rad_m1/implicit_blend_xthin=0 rad_m1/implicit_thin_switch=fse rad_m1/implicit_thin_corr=hcs rad_m1/implicit_thin_recon=jse rad_m1/implicit_thin_bq=lin"
REC="$TJE rad_m1/time_scheme=hesdirk2 rad_m1/time2_vet_sc=predict rad_m1/time2_vstage=true rad_m1/implicit_vimp=false"
RECBE="$TJE rad_m1/time_scheme=be"
# (implicit_vimp=false: kinematic static medium, no gas force; vimp is fatal/inactive there)
run() {  # name, arm keys, extra
  local d=$OUT/$1; mkdir -p $d; cd $d
  $PY $H/fixkeys.py $H/pipe.athinput in.athinput $2 $3 > fixkeys.txt
  echo "$BIN -i in.athinput -d $d $2 $3" > cmd.txt
  OMP_NUM_THREADS=$NT OMP_PROC_BIND=false /usr/bin/time -f "%e s" -o time.txt $BIN -i in.athinput -d $d $2 $3 > run.log 2>&1
  echo $? > rc.txt
}
run REC_c220   "$REC"   "rad_m1/implicit_cfl=220.0" &
run RECbe_c220 "$RECBE" "rad_m1/implicit_cfl=220.0" &
run REC_c2200  "$REC"   "rad_m1/implicit_cfl=2200.0" &
run REC_c22    "$REC"   "rad_m1/implicit_cfl=22.0 time/tlim=1.0" &
wait
for d in $OUT/*/; do echo "$d rc=$(cat $d/rc.txt) $(cat $d/time.txt)"; done
