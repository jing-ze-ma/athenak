#!/bin/bash
# raycheap2-1109 CPU gates (Orion): usage run_gates.sh BIN_NEW BIN_BASE OUTDIR   (MPI CPU builds, PROBLEM none)
# G1 keys absent: new vs base on G_cs_lo_shd_geo2_k8e2 (production K8) and G_cs_lo_shd_geo2_k8e2_ravg_ema16 (1 rank)
#    and on G_cs_lo_shd_ema16_h2 (hesdirk2, no mod keys)
# G2 the flux-mod key acts: G_cs_lo_shd_ema16_h2_mod differs from G_cs_lo_shd_ema16_h2 (new binary)
# G3 restart with the mod key: 20 + 20 vs 40 cycles (1 rank);  G4 mod key 1 vs 4 ranks (cs: 6 MeshBlocks)
BN=$1; BB=$2; O=$3; T=$(cd "$(dirname "$0")" && pwd); mkdir -p $O
run() { local d=$O/$1 np=$2 b=$3; shift 3; mkdir -p $d; (cd $d && mpirun -np $np $b "$@" -d $d > run.log 2>&1; echo "rc=$?" >> run.log); }
run p_new 1 $BN -i $T/G_cs_lo_shd_geo2_k8e2.athinput & run p_base 1 $BB -i $T/G_cs_lo_shd_geo2_k8e2.athinput &
run e_new 1 $BN -i $T/G_cs_lo_shd_geo2_k8e2_ravg_ema16.athinput & run e_base 1 $BB -i $T/G_cs_lo_shd_geo2_k8e2_ravg_ema16.athinput &
run h_new 1 $BN -i $T/G_cs_lo_shd_ema16_h2.athinput & run h_base 1 $BB -i $T/G_cs_lo_shd_ema16_h2.athinput & wait
run m_full 1 $BN -i $T/G_cs_lo_shd_ema16_h2_mod.athinput & run m_p1 1 $BN -i $T/G_cs_lo_shd_ema16_h2_mod.athinput time/nlim=20 &
run m_np4 4 $BN -i $T/G_cs_lo_shd_ema16_h2_mod.athinput & wait
run m_p2 1 $BN -r $(ls $O/m_p1/rst/*.rst | sort | tail -1) time/nlim=40
for d in $O/*; do echo "$(basename $d) $(tail -1 $d/run.log) FATAL=$(grep -a -c FATAL $d/run.log)"; done
C="python3 $T/cmpbin.py"
echo "G1 prod:   $($C $O/p_new $O/p_base | tr '\n' ' ')"
echo "G1 ema16:  $($C $O/e_new $O/e_base | tr '\n' ' ')"
echo "G1 h2:     $($C $O/h_new $O/h_base | tr '\n' ' ')"
echo "G2 mod acts (must DIFFER): $($C $O/m_full $O/h_new | tr '\n' ' ')"
echo "G3 restart: $($C $O/m_p2 $O/m_full | tr '\n' ' ')"
echo "G4 1 vs 4:  $($C $O/m_np4 $O/m_full | tr '\n' ' ')"
