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
