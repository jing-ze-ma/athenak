---
name: ck-nquad2-production
description: User decision 09-25 -- dhj production uses problem/ck_nquad = 2 (two-point Gauss, exact thick limit); nquad 1 (mu=1/1.66) is 10 % low in the diffusion limit
metadata:
  type: project
---
User 09-25: "it matters, no need to test that. make sure nquad=2 works." ck_nquad = 1 (hemispheric mean,
mu = 1/1.66) gives a thick-limit flux 2/D = 1.20 vs the exact 4/3 (tst/test_suite/rad/test_rad_dhj_ck_cpu.py); ck_nquad = 2
(half-range Gauss mu 0.2113/0.7887) is exact. Validation of nq2 with implicit T4 + c2 + every=4: branch ck-nq2,
/viper/ptmp2/jinma/nq2_0925. The ck RCE initial profile uses nq2 (ckrce_nq2.txt, [[sparc-sponge-campaign-0925]]).

**How to apply:** every new dhj input sets problem/ck_nquad = 2; expect roughly 2x the ck cost (measure on GPU).

**Update 09-25 (user):** the ck/Rosseland handover mismatch (two-stream flux 7 % (nq1) / 18 % (nq2) above the table
Rosseland diffusion at 0.02 < w < 0.9; the table kappa_R is not the ck effective mean) is to be removed by option B:
ck two-stream nq2 all the way to the inner wall, no Rosseland blend in x1 -- "go with B if it works"; fallback A =
band-resolved diffusion identical to the two-stream's thick limit. Agent on branch ck-nq2.

**CORRECTION 09-25:** the "7 % / 18 % handover mismatch" was an artifact of the tst check's formula (upper cell's T^3
on a coarse grid), NOT of kappa_R. The table kappa_R matches the exact ck harmonic mean to 0.7-1.3 %; nq2 two-stream /
exact non-grey diffusion = 1.0000-1.0010; handover to the conduction flux ~0.98 (nq2) vs 0.88 (nq1). ck-nq2 f232aecb
(not merged yet): dead-cell fix (CkDeadE, identity rows; both nquad), nq2 speed steps (beam/kappa per angle pair,
lin1p both angles per thread, lP slots once), route B = isotropic_conduction = none + ck_pcut_bar = 1e8. GPU before
speed steps: RT 12.57 (nq1) vs 18.25 (nq2) ms/cycle = 1.45x, total 28.0 vs 33.7 ms/cycle (C32, nx1 128, 2 GPUs).
Pending apudev jobs 11969285/7/8/9 (speed build timings, accuracy + B check).

**ck table extended 09-25 (branch ck-hitemp 8aadfc7b, on ck-nq2; not merged):** Premixed_1x_g8_11_hiT.txt + CE_tables/
FastChem_ck_1x_int_hiT.txt to 10100 K (physics extension: FastChem 3.1.3 composition, line k scaled by neutral carriers,
H bf+ff added, H- rising); bitwise below 6100 K; kappa_R(230 bar) 64->86 (6500 K), 72->231 (8000 K). New key
problem/ck_ce_table; ck_table + ck_ce_table must be IN the input file. Files gitignored: regenerate with
tools/gen_hitemp.py (needs `pip --user pyfastchem`; on orion too). No 12th band (band 10 carries the folded blue tail:
2.8 % of B at 6500 K, 18 % at 10000 K). Merge order: ck-nq2 (after GPU gates) -> ck-hitemp -> dhj-cksph-guard.

**Next (user 09-25):** make implicit nq2 ck (route B, deeper columns) faster AND scale well -- launch a profile-driven agent (branch ck-fast2 on ck-hitemp) once the nq2 GPU jobs 11969285/7/8/9 are analysed and ck-nq2 merged.

**MERGED 09-25:** ck-nq2 d3a04acc, ck-hitemp db8b4764 (pushed). ROUTE B RECOMMENDED for production: <hydro|mhd>
isotropic_conduction = none + problem/ck_pcut_bar = 1e8 + ck_nquad = 2 + T4 + c2 + every 4 (+ hiT tables via ck_table/
ck_ce_table in the input). Accuracy vs tight refs (300-cycle chaotic run, rms eint): b2 1.1e-2, n2 1.0e-2, n1 1.6e-2;
B vs blend 3.5e-3; 0 non-converged. nx1 256 timing job 11969850 pending (analyse: nq2_0925/gpu `python3 ana.py time t256b`).
nx1 256 (MHD C32, 2 GPU): no RT 31.3, n1 76.8, n2 81.3, B 75.3 ms/cycle; B = 0.97x n1, 0 nonconv, 3.0 passes. ck_spherical guard merged d38ace15 + test inputs flipped 9de6060b.

**ck-fast2 STARTED 09-25 (user):** branch ck-fast2, /viper/ptmp2/jinma/ckfast2_0925, target RT <= 0.3x hydro at nx1 256 (today B: 1.0x at 128, 1.4x at 256). Levers in order: (1) consistent deep handover to band-resolved diffusion from the same ck face quantities where all chains are thick, (2) multi-rate cadence, (3) kernel work. USER: KEEP 11 bands x 8 g-points and nquad 2 -- no fewer g-points.
ck-fast2 lever1 09-25: ck_dif_dtau (68a56a44) accurate (flux 1.000000 at handover, eint rms <=6e-8) but only RT -7 %/-2.6 % (nx1 128/256; only 27-36 thick cells/column). Beam chord walk O(N^2) = 52 % of the storing chain kernel at nx1 256 -> lever 3a ck_beam_par. Lever 2 early: d10_e8_x16 RT/h 0.37 at nx1 128, 0 nonconv (accuracy pending).
ck-fast2 lever 2 09-25: e8 (ck_impl_every 8) + x16 (ck_impl_xstep 16) -> RT/h 0.337 at nx1 256 (0.375 at 128), 0 nonconv, accuracy ~1.3x route B's own lever error (eint rms 3.7e-5 at 1-10 bar ... 4.6e-2 top); e16/e32 rejected. Lever 3a ck_beam_par 796a92a5 CPU bitwise, beam 387 -> 172 ms per 8 calls at nx1 256.
ck-fast2 lever 3a 09-25: beam_par (796a92a5+b03f9790) -25 % RT; d10_e8_x16_bt RT/h 0.294 at nx1 256 (total 0.56x b2), 0.36 at 128; BUT e8 leaves 1 stuck column per ~8 calls (9-19 nonconv/12000 s even maxit 12); x16 alone (every 4) 0 nonconv, b2 accuracy, RT/h 0.59. dtau 10 recommended (3 too loose). Fresh start from ckrce: b2 itself 20 nonconv (generic fresh-start issue). Running: guard/e8_x8/e6_x12/maxit20/e4_x16.
ck-hitemp2 merged 09-25 18:00: smooth *_hiT2 tables (1x repo data dir, 10x wasp121_0925/ckdata10, K20 via symlinks); w121/k20 inputs point at _hiT2 (top7/top8 still _hiT). Slope-jump cause was the H bf/ff term added as kH(T)-kH(6100) clipped. Does NOT fix the 10x Newton failure.
no-accuracy-sacrifice rule also applies to ck-fast2 cadence levers (told the agent 21:37).
ck-fast2 MERGED 09-25 ~22:20: recommend ck_dif_dtau=10 + ck_beam_par=true (RT/h 0.79 @128, 0.97 @256); cadence e8/x16 rejected (stuck columns; guard every_thr 0.05 fixes but RT/h 2.1). Not yet in the sponge-run binary.

**User 09-25 ~22:45: 10x ck Newton = option (a): ck_impl_rsec 20 + ck_impl_maxit 12 + ck_impl_tol 1e-7 (user-approved EXCEPTION to the no-accuracy-sacrifice rule, 10x only; measured max residual 1.6e-7).** Applied to wasp121 and kelt20 10x inputs; new binary from rt-integration HEAD in /viper/ptmp2/jinma/bin_0925 (sparc_0925/athena.gpu untouched for the 1x sponge runs).
