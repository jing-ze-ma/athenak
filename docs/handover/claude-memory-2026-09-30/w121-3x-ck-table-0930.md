---
name: w121-3x-ck-table-0930
description: 3x-solar premixed ck table for WASP-121b built 09-30 by ln(k/X) interpolation in log Z between Exo-FMS 1x and 10x; validated on Sonora 2020 at exactly 1 dex; files /viper/ptmp2/jinma/ck3x_0930/ckdata3
metadata:
  type: project
---
No 3x premixed table exists upstream (Exo-FMS: 1x/10x/100x/1000x only; README "more upon request"). Sonora 2020
(Lupu, Freedman, Visscher; Zenodo 10.5281/zenodo.5590997, CC-BY) has the SAME 11 Kataria bands and 4+4 g at
[Fe/H] 0-1.0 incl. 0.5, but T <= 4000 K, P >= 1e-6 bar and different line lists (IR 1.1-1.7x Exo-FMS) - used only as
a validation. Method: interpolate ln(k/X) (per H nucleus; plain ln k is 25-30 % low because mu changes) in log Z;
CE: mu linear in Z, ln VMR linear in log Z. Validation: Sonora 0+1.0 -> 0.5 grey Planck/Rosseland median 1-3 %,
90 % < 8 %; Exo-FMS 1x+100x -> 10x (2 dex) median 8-12 %, 90 % 22 %. Worst: optical bands and 3.5-4.4 um (CO2).
Files ckdata3/ck/Premixed_3x_g8_11_hiT2.txt + CE_tables/FastChem_ck_3x_int_hiT2.txt (hiT ext via gen_hitemp.py,
FastChem +0.47712 dex). EOS: eos_xh 0.7188, eos_yhe 0.2420, problem/met 0.4771. Caveat: e-/H- up to +18 % at
5000-6100 K. Purpose: test whether ~3x metallicity reproduces the observed day-night transport without drag
([[w121-observables-lit-0929]]).

**3x PRODUCTION SUBMITTED 09-30 ~01:15** (w121prod_0929/w3x, jobs 12030240 + afterany 12030241, apu 1 node, TROT 300; binary athena.gpu.w3x = rt-integration 2fd94098; note 1df6f593): EOS met3/dump, RCE IC (T10/T100 = 2627/3831 K, RCB 5.7 bar), grav 1172.79, x1 1.124210e10-1.595088e10, nx1 76 (user 09-30: tool fit; resubmitted 12030356 + 12030357). Smoke clean, dt 16.5 s.

USER 09-30 ~01:30: 3x also handed to Caltech (TASK-2026-09-30-caltech-w121-3x.md 0b893ef6; package fork branch data-w121-3x-pkg, 2.8 MB md5 9dc9472e); viper copy 12030356 stays queued; first to start runs, user cancels the other.

USER 09-30 ~01:45: keep BOTH 3x copies for now (viper 12030356 + Caltech 3640301/02); no cancel instruction to Caltech.

USER 09-30 ~03:10: viper 3x copy CANCELLED; the 3x runs on Caltech only (3640302, started 09-29 17:46 PDT, binary 793e03c3, dt ~12 s, ~24.6 rot/h, ETA rot 300 = 09-30 06:00 PDT / 15:00 CEST; dir /resnick/groups/carnegie_poc/jingze/w121prod_0930/).

**09-30 ~04:30 Caltech 3x dt decline 12 -> 7.5 s (rot 18-30), radial CFL from supersonic night-side downflows (64 columns, to 1e-5 bar); hypothesis: ck solver limit cycle with 3x/1x keys (dtmax 0.5, maxit 16, no rsec) vs 10x keys (0.25, 24, rsec 20, tol 1e-7). Caltech A/B from rot-36 rst (3649979 3x keys / 3649980 10x keys). USER: wait for it; viper repeat only if ambiguous.

USER 09-30 ~05:25: viper 3x with ONLY ck_impl_dtmax 0.25 + maxit 24 (w3xk, fresh start; tol 1e-8, no rsec): jobs 12034146 (+12034147); input w121prod_3x_ck10.athinput; binary cba4e797; compare with Caltech 3x (3x keys) on dt decline / night-side downflows / top-cell cycle.
w3xk final jobs: 12034178 (+12034179)

USER 09-30 ~08:00: w3xk (12034178/79) on HOLD until the lhllc radial-face fix is ready (then add the chosen key to w121prod_3x_ck10.athinput + rebuilt binary, resubmit or release).

USER 09-30 ~09:35: w3xk RELEASED with the radial-fix binary (lhllc-x1-phi 36dd494d, md5 30d02515) + lhllc_x1_phi_min = 1 in its input; smoke clean. So w3xk = 3x with dtmax 0.25 + maxit 24 AND the radial fix.

USER 09-30 ~09:55: Caltech A/B: no ck limit cycle with 3x keys (dt 7.84 vs 7.80 s) -> dt decline not a solver effect. w3xk input RESET to the original 3x solver keys: w3xk = 3x + radial fix only. Fixed set w1xf/w3xk/w10xf vs pre-fix w1x/Caltech 3x/w10x (watch the 3x dt decline).
