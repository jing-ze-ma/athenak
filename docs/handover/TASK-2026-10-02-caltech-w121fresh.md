# TASK for Caltech: WASP-121b fresh starts 1x / 3x / 10x on face 6 + wall-flux fix (duplicates of viper's queued chains)

User decision 10-02 07:20: run the three viper WASP-121b fresh-start arms on Caltech as well; for each arm whichever site
STARTS first is kept and the other is cancelled (viper does the cancelling on its side; push a NOTE the moment one of
yours starts). Viper chains (pending): f1x6w 12059646->47, f3x6w 12059648->49, f10x6w 12059650->51->52.

## 1. Code
- Build the dhj problem (deep_hot_jupiter_rt, your usual WASP-121b CUDA build recipe) from **rt-integration at or after
  3f3b0414** (this commit is fine). It contains, all as fresh-run defaults (logged at startup, please confirm):
  ck_sph_face = 6 (monotone ladder; + companions ck_sph_top 1, ck_vef_every 150, ck_impl_rowfb, ck_impl_rsec 4),
  <problem>/ck_wall_flux_exact = true (faces 5/6 injected 7.2 % too much internal flux before), the AM fix
  (coriolis_am / am_transport when rotating), radial top sponge, implicit_face_weight distance (M1 only, unused here).
  Your EXISTING WASP-121b runs predate face 6 and the wall fix: do not restart them with this binary for this task.

## 2. Inputs (docs/handover/w121fresh_1002_bundle/)
- w121fresh_{1x,3x,10x}_f6.athinput = viper's production inputs with only the data paths replaced by placeholders:
  @BUNDLE@ (the bundle dir: IC files ic_w121_*_f5.txt, MD5_IC), @CK1X@ / @CK3X@ / @CK10X@ = the Exo-FMS ck data dirs for
  1x / 3x (met 0.4771) / 10x (met 1.0) with ck/Premixed_*_g8_11_hiT2.txt inside. The tables are NOT in git (14 MB each):
  use your local copies and check them against MD5_CK (paths CK1X/..., CK3X/..., CK10X/...). If any md5 differs, STOP
  and report (do not substitute another table).
- The ICs are the exact-reference (p-z formal solution) RCE profiles from viper (w121fresh_1001/ic/rce_f5.py), 90-160 K
  warmer at 1e-3..0.1 bar than the old ICs. Grid/physics otherwise = your earlier w121prod inputs (cubed sphere C32,
  lhllc, etc.); time/tlim = 300 rotations as in the input.

## 3. Gates
- Startup log shows ck_sph_face = 6, ck_wall_flux_exact = 1 and the companions; a short CPU vs GPU check if your recipe
  has one; then a smoke (fresh start ~5 min + a restart leg): viper's smoke dt medians were 15.2 s (1x), 12.9 s (3x),
  12.6 s (10x), minimum ~8-13 s; 0 FATAL/nan. The 10x arm is the multi-rotation check of face 6 at high metallicity:
  per-rotation dt median >= 5 s, KE < 1e34, |Etot/Etot0-1| < 1e-2 (viper check_xl.py logic).
- Your GPU count per node differs from viper's 2 MI300A, so the MeshBlock layout differs: runs are statistically, not
  bitwise, comparable. Fine for this relaxation run.

## 4. Production
- Chained 24 h links per arm to 300 rotations with a link sanity check (refuse on NaN/FATAL, no newer rst, or median dt
  < 2 s), restarts with every <outputN>/last_time reset to floor(t/dt)*dt. Outputs as in the input.
- Report: NOTE-2026-10-0x-caltech-w121fresh.md at submit (job ids, layout, s/cycle) and at START of each arm
  (so viper cancels its copy); then every ~50 rotations (dt, OLR/absorbed, deep T at 10/100 bar, any FATAL).
