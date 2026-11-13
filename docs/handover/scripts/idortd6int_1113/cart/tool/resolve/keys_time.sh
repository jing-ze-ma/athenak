# TIME-ONLY candidate (physflux-1111 86714ac4): production deck keys + k1_redist 0 (as BASE/CAND of the 11-11 gate)
# + cmem + cap 1 + x2f15 + physflux + physflux_mom fa; NO space keys (no plm-on-h, sfgate, gjump, corr_plm, (a), (b), sw5)
PT="rad_m1/implicit_hr_cmem=true rad_m1/implicit_hr_cmem_cap=1.0 rad_m1/time2_vet_x2_fields=15 rad_m1/implicit_physflux=true rad_m1/implicit_physflux_mom=fa"
KT="rad_m1/time2_k1_redist=0 $PT"
# curvilinear: the 11-11 KS minus its space keys (recon plm/q/lim, thin_recon dc, gface jump, corr_plm, h0/h1 sw5, hlim, qe_tau)
KTS="rad_m1/vet_col_ncore=64 rad_m1/vet_gd_upwind=shellmean rad_m1/vet_gd_upwind_nbin=8 rad_m1/vet_gd_rotate_every=0 rad_m1/vet_gd_rotate_dt=0 rad_m1/vet_rotate_avg=1 rad_m1/vet_rotate_avg_ema=0 rad_m1/vet_gd_upwind_wall=true rad_m1/implicit_thin_corr=hc $PT"
# idort-x1-1112: X1 additions (combo-c3-1112 559003d18) on top of KT / KTS (time2_vet_x2_fields=15 already in PT)
X1="rad_m1/implicit_physflux_zold=fc rad_m1/implicit_hr_cmem_rst=true rad_m1/implicit_hr_cmem_fainit=true rad_m1/time2_stage_local_be=true rad_m1/time2_face_fct=true rad_m1/time2_k1_redist=0"
# + rad_m1/vet_gd_rotate_xfade=true only where vet_gd rotation is on the physical clock after all keys (tool/xfade.py decides)
XF="rad_m1/vet_gd_rotate_xfade=true"
# idort-x1-1112 FC (USER 11-12: the suite runs the FINAL candidate FC). Two candidates, chosen by CAND (env, else $W/cand.txt):
#   cmem = cmemimpl-1112 aea1cae35 (cmem v2; v1 7bc19d109 defective), key implicit_hr_cmem_implicit
#   fmm  = facemmat-trans-1112 c6e42d5f1, key implicit_face_mmat (the two keys FATAL together in the cmem code)
# FC = X1 + FCA on every deck + FCC on sp / cs decks (all read only when named). Binaries: $W/bin_{mpi,omp}.<CAND>.txt
CAND=${CAND:-$(cat /viper/ptmp2/jinma/d6int_1113/orion_task/cand.txt 2>/dev/null)}
case "$CAND" in cmem) FCA="rad_m1/implicit_hr_cmem_implicit=true rad_m1/implicit_cs2_wbsolve=true";;
                fmm)  FCA="rad_m1/implicit_face_mmat=true rad_m1/implicit_cs2_wbsolve=true";;
                fin)  # USER/coordinator 11-12 FINAL (RSG verdict), combo-all-1112 421a3d024: X1b = TIME + stage_local_be + face_fct +
                      # k1_redist 0 (+ xfade where vet_gd rotates; x2_fields 15 is in PT); NO physflux_zold/cmem_rst/cmem_fainit, NO cmem_implicit
                      X1="rad_m1/time2_stage_local_be=true rad_m1/time2_face_fct=true rad_m1/time2_k1_redist=0"
                      FCA="rad_m1/implicit_face_mmat=true rad_m1/implicit_cs2_wbsolve=true rad_m1/implicit_am_cs_fold=applied";;
                d6)   # D6 (d6gates 11-13): fin minus implicit_am_cs_fold, cmem/physflux off (PT sets them on), blend_dtc 1e-6
                      X1="rad_m1/time2_stage_local_be=true rad_m1/time2_face_fct=true rad_m1/time2_k1_redist=0"
                      FCA="rad_m1/implicit_face_mmat=true rad_m1/implicit_cs2_wbsolve=true rad_m1/implicit_hr_cmem=false rad_m1/implicit_physflux=false rad_m1/implicit_blend_dtc=1e-6";;
                *) echo "keys_time.sh: CAND='$CAND' is not cmem|fmm" >&2; FCA="CAND_UNSET";; esac
FCC="rad_m1/implicit_cs_e0_metric=true rad_m1/vet_col_lat_posfix=lag rad_m1/implicit_comoving_lat=true rad_m1/implicit_cs_closure_metric=true"
# arm p = CAND (cand.txt); optional second arm q = always fmm (keys + binary), for both candidates side by side: ARMS="p q", cs ARM=q
# cand fin: FCC (sp/cs) = vet_col_lat_posfix=lag + implicit_cs_e0_metric only (no comoving_lat / closure_metric)
{ [ "$CAND" = fin ] || [ "$CAND" = d6 ]; } && FCC="rad_m1/implicit_cs_e0_metric=true rad_m1/vet_col_lat_posfix=lag"
