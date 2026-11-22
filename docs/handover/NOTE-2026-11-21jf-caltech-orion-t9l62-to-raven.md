# NOTE 2026-11-21jf Caltech -> orion: USER — continue t9l62 on RAVEN (from the Caltech final rst)

USER (in our session, 17:2x PST): "continue t9l62 on raven". Everything needed:

- Restart file: Caltech hpc /resnick/groups/carnegie_poc/jingze/prod_fresh_1121b/t9l62/rst/t9l62.00037.rst
  (27,988,115,695 bytes, md5 6c839345c394ecb99234085b3725f591, cycle 30085, t = 0.587 tlim; written by the B_p1 wall-time stop, rc 0).
- Deck: J4_t9l62.athinput (copy: caltech_results/PROD1121b/t9l62_J4_deck.athinput); restart reads its keys from the rst.
- Code: caltech-prod-bf4f917d4 (B_p1). Build for A100 on Raven (CUDA 12.x >= 12.9 preferred, never 12.2 per the H200 trap; smoke-gate on Raven as you do).
- Command-line keys (what Caltech ran since 16:41): rad_m1/clean_conv_max=1.0e-4 rad_m1/clean_face_opac=harmonic rad_m1/clean_stag_opac_update=true
  rad_m1/clean_stag_flim=false rad_m1/clean_stag_fwb=true rad_m1/clean_stag_fcone=true rad_m1/clean_stag_fwb_stage=true rad_m1/clean_stag_thinf_held=false
  rad_m1/clean_gas_newton_switch=4 rad_m1/clean_gas_newton_switch_min=6 time/cfl_number=0.6 ; env VGD_XFE=0. No damp, no xpat. tlim 2.682e6 s.
  (cfl 0.6: at 0.9 it had a stage fallback within 100 cycles at 16:39; the 0.6 chunk had CW 1, FATAL 0, retry 0.)
- Layout on Caltech: 1 node x 4 H200, NR=4, 8 MeshBlocks 512x128x128 (~80 GB/GPU peak on H200): on A100-40GB this needs more GPUs (e.g. 2 nodes x 4 A100 = 8 ranks, one MeshBlock each); A100-80GB 4 GPUs may fit.
- TRANSFER: Caltech cannot reach raven.mpcdf.mpg.de (ssh port 22 times out from our side). Please say how you want the 28 GB rst moved (you pull from Caltech, Globus, or via an intermediate host) and we act at once. The rst stays in place; it is hard-linked so the pruner cannot delete it.
Caltech stops t9l62 here (holders end 17:46); it is not relaunched on Caltech while it runs on Raven.

-- Caltech
