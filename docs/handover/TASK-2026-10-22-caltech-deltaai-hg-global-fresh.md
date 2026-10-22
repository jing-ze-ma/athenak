# TASK 2026-10-22 caltech -> deltaai (USER): HG GLOBAL production, FRESH start on DeltaAI with the passing collapse fix

Why: Caltech's HG global (and RSG E) harmonic-only runs hit a dt collapse (HG at 18.8 d, dt 125 -> 1e-4 s). Physical photospheric
clumping (RT-type) + solver rejections. The fix that passed on Caltech (RSG E arms E / E2, NOTEs 01599ea3, f0bd7ea1):
  <rad_m1>/implicit_face_opac = harmonic, time2_vet_col = extrap2, time2_accept_fallback = true, time2_posfix = 4
USER: run HG global FRESH (t = 0) on DeltaAI with that fix, racing Caltech's queue (Caltech continues HG from its 15 d rst; HG 3-node
7-d allocation est. tonight; debug-QOS links running). Report the job ids, start time and s/cycle in a NOTE; USER decides which copy continues.

Code: fork fin51-r-fix = ff1f52bc (Caltech H200 binary md5 cd998ac0, nofma: -DCMAKE_CXX_FLAGS=-ffp-contract=off, PROBLEM=he_star_m1).
Build INCREMENTALLY if you have a configured he_star_m1 tree (USER rule), on GH200 (Kokkos_ARCH_HOPPER90).
Deck: docs/handover/hgglob-deltaai-1022/hgglob_fresh_clip_harm.athinput (= Caltech deck md5 above, file paths -> @FILES@):
  replace @FILES@ with your checkout's docs/handover/hg145-prod-1018 (ic_hegiant_14p5M_mk.txt f9678c2a, rosseland_mesakap_afe2dd8da2bd.txt
  e2868b75, planck_ext2_..._416a61196c64.txt 23ee543f -- all already on bsg-files-1009).
Command-line keys (the Caltech production set):
  rad_m1/time2_vet_col=extrap2 rad_m1/time2_accept_fallback=true rad_m1/time2_posfix=4 output2/dt=28800 output3/dt=28800 time/tlim=1.0e9
Layout as Caltech: 24 MeshBlocks 704x64x64, 3 nodes x 4 GPUs (12 ranks, 2 MB per GPU), or what fits your queue best; longest walltime links.
Do NOT use <hydro>/am_vector (posfix not wired for it; TASK 92c087ad to viper). Pages: Caltech HG URL LNtjtaL34oUhgcCRNp5V2M is Caltech's run.
