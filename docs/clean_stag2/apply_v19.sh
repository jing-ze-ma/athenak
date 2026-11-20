#!/bin/bash
# apply_v19.sh (INTEG): the USER's v19 key decision (2026-11-19 night, conditional on the evidence, met):
#   rad_m1/time2_face_fct_sweep=true (USER approved; CPU bench gate/tab_swS_vs_tolD.md: 66 rows =, FCT rounds ~ -50 %; MI300A HG -10 %) for HG and RSG (not R136) in the LATEST keys files + integ k1 args. Backups *.bak_v18.
#   (The A100 RSG -19 % / tab_tolK_vs_tolD.md evidence was for clean_conv_kappa, which is NOT in v19.)
#   Run ONLY when v19 is cut; the Raven speedhead keys (A100HR sh_*) are updated at the same time, not before.
S=/orion/ptmp/jinma/clean_1115/stag2
# USER decision (coordinator 2026-11-20): v19 = sweep ONLY. clean_conv_kappa=1.0 is OUT of v19 (A100: kappa without the sweep
# cascades, 7b3003931 RSG dt 1000 -> 62.5, 14 CW; with the sweep clean over 75 cycles); revisit after the DeltaAI kappa+sweep long run.
# clean_conv_kappa=1.0 must never be used without time2_face_fct_sweep=true.
# SWEEP=0: without the sweep (coordinator recommends dropping it: MI300A +3-5 %, A100 RSG +1.4 %, no robustness gain; USER decides)
KV="rad_m1/time2_face_fct_sweep=true"; [ "${SWEEP:-0}" = 0 ] && KV=""   # USER 2026-11-20: sweep DROPPED from v19 (default SWEEP=0)
# USER 2026-11-20: rad_m1/clean_rot_fuse=true (HG only; RAYS 3f8e5e37f, default off): A100 HG 12r 39.48 GB, rt bitwise,
# hst = off, HG 24r -2.5 % mid. Not RSG / R136. The HG 3-node variant = hg_latest + time/dt_retry_host=m1 inherits it.
KH="rad_m1/clean_rot_fuse=true"
# INTEG3 2026-11-20 (USER "merge and make some defaults"): on the R1 head (GEOM b05186c50 = INTEG 458df0a0a) clean_cs_lz_abs=2 and
# clean_stag_csvec_radial are the cubed sphere's only form (old values FATAL): NO lz / csvec keys in any keys file.
# clean_conv_kappa stays at the code default 0.1 (no key). Sweep off (SWEEP=0). HG clean_rot_fuse=true (KH).
for f in $S/hg_latest.txt $S/rsg_latest.txt; do cp $f $f.bak_v18; K="$KV"; case $f in *hg_latest*) K="$KV $KH";; esac; for kv in $K; do k=${kv%%=*}; grep -v "^$k=" $f > $f.t && mv $f.t $f; echo "$kv" >> $f; done; done
for f in $S/integ/keys/k1_hg.args $S/integ/keys/k1_rsg.args; do cp $f $f.bak_v18; K="$KV"; case $f in *k1_hg*) K="$KV $KH";; esac; for kv in $K; do sed -i "s|[ ]*${kv%%=*}=[^ ]*||g" $f; done; [ -n "$K" ] && echo " $K" >> $f; done
grep -H "clean_conv_kappa\|face_fct_sweep\|clean_rot_fuse" $S/hg_latest.txt $S/rsg_latest.txt $S/r136_latest.txt; grep -o "clean_conv_kappa=[0-9.]*\|face_fct_sweep=[a-z]*\|clean_rot_fuse=[a-z]*" $S/integ/keys/k1_*.args
