# beams-1110: the Cartesian beam suite of hr2prod-1110 (from Raven /raven/ptmp/jinma/hr2prod_1110/beams and
# /raven/ptmp/jinma/order_gpu_1110/beams/ana, copied 11-11)
- decks/: cb_FIN3_r4, ws66, shd_FIN3_r12g2, hcyl_FIN3_g8x32 (blend_1022 decks). Every arm adds rad_m1/time2_k1_redist=0.
- run_beams.sh: how the arms were run on Raven (A100); its key() table maps arm names to command-line keys. For the
  Cartesian cand / w1 / wc5 add the TASK keys (setkeys.py for keys absent from a deck). Paths inside are Raven paths.
- ana/score.py ROOT [deck ...]: scores ROOT/<deck>/<arm>/rep<r>/ (run.log, rc.txt, time.txt, in.athinput, bin/) ->
  cb E_L1, ws66 E_L1, shd E_L1, hcyl E_out / rF. Base (HR2ND_PROD 5): 0.5700 / 1.1771 / 0.0131 / 0.0138, 1.0081.
  Exact references: ana/*_exact_*.npz.
