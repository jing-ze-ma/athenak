#!/bin/bash
# run_m35.sh A|B GPUS: orion NOTE 23zz9 step 2: IA12.2 (athena_ia122_465aff485) m35mms 300-cycle restart arms from prod_ia11_1123/m35mms/rst_keep/m35mms.00014.rst (cycle 33724, read only)
# -i ia122_decks/notile/IA122_m35mms.athinput (production m35mms layout = block -> notile) names every IA12.2 key; B adds vet_gd_wall_sym + vet_gd_theta_seam_sym; no rst/bin output
J=/resnick/groups/carnegie_poc/jingze; A=$1; GP=$2; D=$J/ia122_1123/m35_$A; BIN=$J/wt_orion1121/bin/athena_ia122_465aff485; RST=$J/prod_ia11_1123/m35mms/rst_keep/m35mms.00014.rst; DK=$J/ia122_1123/decks/notile/IA122_m35mms.athinput; MW=$J/orion_runs_1121/ctl/memwrap.sh
EX="time/nlim=34024 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30"; [ $A = B ] && EX="$EX rad_m1/vet_gd_wall_sym=true rad_m1/vet_gd_theta_seam_sym=true"
P=$J/prod_fresh_1121b; while IFS= read -r l; do export "$l" 2>/dev/null; done < $P/i1_m35mms_f1k.env; export VGD_XFE=0
mkdir -p $D; cd $D; echo "$(date '+%F %T') start m35_$A BIN md5=$(cut -c1-32 $BIN.md5) deck md5=$(md5sum < $DK | cut -c1-32) GPUs=$GP RST=$RST EX: $EX" >> $D/run.log
OMPI_MCA_io=romio321 srun --overlap --jobid=7333072 -w hpc-sm-01-17 -N 1 -n 4 --mpi=pmix -t 120 -J m35$A \
  env CUDA_VISIBLE_DEVICES=$GP KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $MW $D $BIN -r $RST -i $DK -d $D $EX < /dev/null > $D/log.0 2>&1
echo "$(date '+%F %T') rc=$? last: $(grep -E '^elapsed=.*cycle=' $D/log.0 | tail -1)" >> $D/run.log
