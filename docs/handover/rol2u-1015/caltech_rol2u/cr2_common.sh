# cr2_common.sh (CR2 caltech_rol2u_1015): sourced by seed_cr2.sh / gate_cr2.sh / win_cr2.sh. Caltech paths; edit here only.
G=/resnick/groups/carnegie_poc/jingze
W=${CR2_W:-$G/rol2u_1015}                       # CR2 work dir (NEVER a production dir)
X=${CR2_X:-$W/caltech_rol2u}                    # this bundle (cr2ana.py, cmprst.py, bin_convert.py)
SCR=$G/e1s_1014/scripts                         # rst_pick.py, cmplast.py, bin_convert.py (e1caltech-1013 set)
BIND=/resnick/home/jingze/ATHENAK/builds/bin
BQ=$BIND/athena_gpu_he_star_m1_1bd0c2ce_nofma   # fin51-r @ 1bd0c2ce (build per TASK section 1)
BQ_MD5=${BQ_MD5:?set BQ_MD5 = md5 of the 1bd0c2ce nofma build (TASK section 1)}
ROL="rad_m1/implicit_tol_rel=0.1 rad_m1/implicit_tol_rel_max=1e-4 rad_m1/implicit_opac_lag=pred rad_m1/implicit_opac_lag_passes=2 rad_m1/implicit_gas_newton_u4=true"
DR="time/dt_retry=true"                         # E1 deck has it already (harmless); HG FIN4 rst lacks it (absent key: added on restart)
case ${CASE:?CASE=hg|e1} in
  hg) PROD=$G/hegiant_1014/wide4/run; PFX=hegiant
      BP=$BIND/athena_gpu_he_star_m1_d505d21f_nofma; BP_MD5=1440c05320a7bef88c4ad960019b9840
      NGPU=4; LAUNCH=${LAUNCH_HG:-$G/e1s_1014/launch4.sh}   # = the hgF4w launch (1 rank / H200, 1 MeshBlock 897x64x64 each)
      PKEYS="time/tlim=1e9"; NG=30; WIN=2.5e5;;
  e1) PROD=$G/e1s_1014/prod; PFX=agcar3d
      BP=$BIND/athena_gpu_he_star_m1_5b7e1bc5_nofma; BP_MD5=6252e58e5192ce65ff1f1f627e1ce90c
      NGPU=2; LAUNCH=${LAUNCH_E1:-$G/e1s_1014/launch2.sh}   # = agcE1s launch (2 ranks x 2 MeshBlocks 536x64x64)
      PKEYS="time/tlim=1e9"; NG=20; WIN=2.0e6;;
esac
D=$W/$CASE; SEED=$D/seed
export UCX_TLS=^xpmem OMP_NUM_THREADS=1
chkmd5() { local m=$(md5sum "$1" | cut -d' ' -f1); [ "$m" = "$2" ] || { echo "md5 $1 = $m != $2"; exit 1; }; }
seedrst() { ls $SEED/rst/*.rst 2>/dev/null | head -1; }
seedkeys() { python3 $SCR/rst_pick.py $SEED/rst | tail -1; }        # outputN/last_time = production output cadence
seedinfo() { python3 $SCR/rst_pick.py $SEED/rst | head -1; }        # "<path> <t> <ncycle> <tlim>"
lastrst() { ls $1/rst/*.rst 2>/dev/null | sort | tail -1; }
parend() { grep -abo '<par_end>' "$1" | head -1 | cut -d: -f1; }   # must stay < 45056 (parameter_input.cpp 11 x 4096-byte reads)
