---
name: ck-spherical-closure-overemits-0930
description: 09-30 emission post-processor finding - the dhj ck_spherical thermal two-stream closure (S=(I+ + I-)/2 continuous at faces) over-emits OLR ~25 % (1x) / ~15 % (10x) for the given T; unverified by a unit test; user decision on a fix
metadata:
  type: project
---
Emission agent (synth2_emis, 09-30 ~23:50): ray RT through the SAME columns with the GCM's own 11-band table gives
L = 0.799 x rt_Fb OLR (1x; 0.872 at 10x), uniform over bands and day/night; re-implementing the GCM ck_spherical
closure (src/utils/two_stream_rt.hpp:988-1013, face conditions S=(I+ + I-)/2 and A D continuous) reproduces rt_Fb to
0.960 (1x) / 1.023 (10x). Mechanism (agent): S continuity across the area jump of the transparent upper shell drives
the inward stream negative, so the photosphere radiates by 2x/(1+x), x = A_top/A(R_p) = 1.715 -> 1.263 predicted,
1.25 measured. Implication: the GCM planet is too cold for its absorbed energy (~5-6 % in T); explains part of the
dayside deficits vs data (NRS2 new 3830 vs obs 4924 ppm, SOSS o1 681 vs 1156). ck_spherical defaults TRUE since
2026-09-22 on sp/cs meshes -> all WASP-121b runs have it. NOT yet checked by a clean unit test (isothermal sphere with
transparent shell: L must equal 4 pi R_ph^2 sigma T^4). A fix changes production thermal structure: user decision.
Related: [[w121-1x-rot300-0929]], [[dhj-rsolver-lhllc-lhlld]].

**CONFIRMED 10-01 ~00:20 (unit test in the REAL kernel, /viper/ptmp2/jinma/cksph_test_0930/RESULTS.md):** sharp photosphere
under a transparent shell emits 2x/(1+x) sigma T^4 to 1e-4 (x = 1.2/1.7/3 -> 1.091/1.259/1.500); plane-parallel emits x.
Static isothermal W121 1x grid: L_code/L_exact 1.224 (2000 K) / 1.248 (3000 K); old dhj grid 1.236/1.338. Cause: face passes
S and A(u-d) = Eddington closure K=J/3, J never dilutes. Prototype problem/ck_sph_dilute (d_below = d_above; A u conserved
in a transparent shell) on branch ck-sph-closure-fix (7314aad3 hooks, e7c50db5 fix, 9ad3f06b tests; not pushed): sharp
tests 1.000, isothermal 0.97/0.98 (W121) - remaining few % = radial column cannot carry grazing chords; deep diffusive flux
0.2-0.3 % low. Only tm kernel + ck_implicit; REFUSES ck_impl_lin/jac_lin (production uses them) -> port needed + GPU A/B.
Thin-shell heating error 8-15 % remains (variable-Eddington factor would fix). Hot top-ghost boundary value flagged.
