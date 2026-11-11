#!/bin/bash
# run_beams.sh BIN GPU "deck:arm ..." : Cartesian beam decks (blend_1022/p2) on one A100, sequential; k1_redist 0 on all arms.
BIN=$1; G=$2; SETS=$3
export CUDA_VISIBLE_DEVICES=$G OMP_NUM_THREADS=1
for v in $(compgen -e | grep -E '^(PMIX|PMI_|OMPI_MCA_ess|OMPI_MCA_pmix|PRTE|SLURM)'); do unset $v; done
DK=/raven/ptmp/jinma/hr2prod_1110/beams/decks
O=/raven/ptmp/jinma/hr2prod_1110/beams; SK=$O/setkeys.py
K0="rad_m1/time2_k1_redist=0"
NOCJ="rad_m1/implicit_thin_corr=none rad_m1/implicit_thin_recon=dc"
PLM="$NOCJ rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_lim=smooth"
key() { case $1 in
  base) echo "$K0";;
  gmean) echo "$K0 rad_m1/implicit_thin_gface=mean";;
  gjump) echo "$K0 rad_m1/implicit_thin_gface=jump";;
  plm) echo "$K0 $PLM";;
  plmgjump) echo "$K0 $PLM rad_m1/implicit_thin_gface=jump";;
  plmgmean) echo "$K0 $PLM rad_m1/implicit_thin_gface=mean";;
  nofse) echo "$K0 rad_m1/implicit_thin_switch=none";;
  cm) echo "$K0 rad_m1/implicit_hr_cmem=true";;
  cmc4) echo "$K0 rad_m1/implicit_hr_cmem=true rad_m1/implicit_cfl=4";;
  cmc16) echo "$K0 rad_m1/implicit_hr_cmem=true rad_m1/implicit_cfl=16";;
  cmc64) echo "$K0 rad_m1/implicit_hr_cmem=true rad_m1/implicit_cfl=64";;
  bc16) echo "$K0 rad_m1/implicit_cfl=16";;
  cmnofse) echo "$K0 rad_m1/implicit_hr_cmem=true rad_m1/implicit_thin_switch=none";;
  cmnocj) echo "$K0 rad_m1/implicit_hr_cmem=true rad_m1/implicit_thin_corr=none rad_m1/implicit_thin_recon=dc";;
  cmgjump) echo "$K0 rad_m1/implicit_hr_cmem=true rad_m1/implicit_thin_gface=jump";;
  dtc) echo "$K0 rad_m1/implicit_blend_dtc=1.0e-6";;
  gjumpdtc) echo "$K0 rad_m1/implicit_thin_gface=jump rad_m1/implicit_blend_dtc=1.0e-6";;
  cen) echo "$K0 rad_m1/implicit_flux=central rad_m1/implicit_flux_beam=closure rad_m1/implicit_thin_switch=none rad_m1/implicit_thin_corr=none rad_m1/implicit_flux_faces=x1 rad_m1/implicit_blend=tau_f";;
esac; }
for x in $SETS; do d=${x%%:*}; a=${x#*:}; D=$O/$d/$a/rep1; mkdir -p $D; cp $DK/$d.athinput $D/in.athinput; python3 $SK $D/in.athinput $(key $a)
  t0=$(date +%s); $BIN -i $D/in.athinput -d $D > $D/run.log 2>&1; echo "rc=$?" > $D/rc.txt
  echo "$(( $(date +%s)-t0 )) s wall" > $D/time.txt; echo "$(date -u +%FT%TZ) beam $d $a $(cat $D/rc.txt) $(cat $D/time.txt)" >> $O/verdicts.txt
done
