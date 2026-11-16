#!/bin/bash
# wb.sh BIN: well-balancing level, max v/c_s over t <= 4 of the rhse envelope at rest with
# force_reference = clean and the hydro static well-balancing (~1 min on 16 cores)
B=$(readlink -f $1); R=$(cd "$(dirname "$0")" && pwd)
export OMP_PROC_BIND=false OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=4
WB="rad_m1/force_reference=clean hydro/wellbalance_static=true hydro/wellbalance_static_reconst=true time/tlim=4"
for r in fv1 fv lld cart; do
  d=$R/wb/runs/$r; mkdir -p $d
  if [ $r = cart ]; then (cd $d && $B -i $R/wb/rhse_cart.athinput $WB time/cfl_number=0.3 > log.txt 2>&1) &
  else (cd $d && $B -i $R/wb/rhse_sp.athinput $WB rad_m1/clean_rays=$r > log.txt 2>&1) & fi
done; wait
for r in fv1 fv lld cart; do d=$R/wb/runs/$r
  echo "wb $r max v/c_s $(awk 'NR>2{print $5}' $d/*user.hst | sort -g | tail -1) CLEAN WARNING $(grep -c 'CLEAN WARNING' $d/log.txt)"; done
