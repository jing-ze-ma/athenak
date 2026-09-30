---
name: lhllc-radial-oddeven-0930
description: lhllc (dhj production solver) leaves a deep 2-cell RADIAL odd-even velocity mode undamped (low-Mach phi ~3e-3); plain hllc removes it within 0.5 rot; mode regrows at nx1 256; deep v_r / Mdot / deep enthalpy flux of lhllc runs untrustworthy
metadata:
  type: project
---
Test 09-30 (/viper/ptmp2/jinma/w121prod_0929/oddeven_0930, oe_{L,H,R}.txt; 1x rot 300 restart, same binary cba4e797):
odd-even share (1 = pure 2-cell mode) at >100 / 30-100 / 1-10 / 0.1-1 bar: L lhllc 0.96/0.93/0.91/0.91 unchanged over
5 rot; H hllc 0.10/0.04/0.08/0.21 by 0.5 rot; R (nx1 256, lhllc) regrows 0.21 -> 0.87 at >100 bar within 0.76 rot.
Per-shell Mdot sign alternation L 1.00 everywhere, H 0-0.25. Radial KE 2.0e32 (L) vs 4.9e30 (H): ~97 % of L's radial KE
is the mode. H: deep horizontal v rms 7-11 %% lower (hllc low-Mach dissipation), jet 0.3 %%, dt 13.8 vs 11.9 s.
Mechanism: lhllc_hyd.hpp:120-123 phi = chi(2-chi), chi = max|v_n|/max c ~1e-3 radially -> the velocity-jump
dissipation in p* (:154-156) is ~250x weaker; mass flux keeps it -> undamped velocity checkerboard (pressure-velocity
decoupling of low-Mach fixes). Consequences: deep v_r, per-shell Mdot and deep resolved enthalpy flux of lhllc runs are
this mode ([[w121-1x-rot300-0929]] deepmix baseline, 10x deep cooling); lhlld (MHD) likely the same (not tested).
Proposed fix (not done): full velocity-jump term or a phi floor on x1 (radial) faces only. Relates to user rule
[[dhj-rsolver-lhllc-lhlld]].

**FIX 09-30 ~08:45 (branch lhllc-x1-phi fbe533f3, local; x1phi_0930/RESULTS.md):** keys <hydro>/lhllc_x1_phi_min, <mhd>/lhlld_x1_phi_min (phi = max(phi, phi_min) on x1 faces; 1 = full HLLC radially); key off bitwise. 1x 5 rot: V1 (1) removes the mode like hllc (oe .04-.20, radial KE 8.8e30 vs 2.0e32) while deep horizontal v rises +1..+6 % (still rising; opposite of hllc -9 %), jet -0.6 %, dt same; V2a/V2b (0.05/0.1) leave a residual mode at 0.1-10 bar. Preliminary: V1. Pending 10x + lhlld smoke. Side issue: xL0 final tlim restart truncated (98 vs 141 MB).

USER 09-30 ~09:00: once the 10x check + lhlld smoke pass, make V1 the DEFAULT (lhllc_x1_phi_min = lhlld_x1_phi_min = 1; 0 = old); merge after combined gate; then swap the new binary into the held w3xk and release.

**DEFAULTS AUDIT 09-30 (docs/handover/AUDIT-2026-09-30-default-candidates.md, uncommitted):** A = ck_impl_maxit 8->24 (ceiling; non-converged calls currently stop on unchecked iterates), ck_impl_conserve 0->1, bbot default 0 (cosmetic). B = lhllc/lhlld_x1_phi_min (GLOBAL default would hit He box/He presn where x1 is the convective direction -> dhj-only default or He box check first), vet_col_source relaxed (ONLY on he-presn-m1, not rt-integration; needs merge + T-S4/Milne), ck_impl_dtmax 0.25 (pending A/B; 1x evidence against). Mismatch: sponge_bottom default true vs production false. Already default: cs_seam_flux positive, wall_closed, rotpot fix, flux_hst_rkavg, opac Newton + guard, M1 2nd-order set, ck T4 set.

USER 09-30 ~09:05: radial fix = dhj-ONLY default (pgen sets 1 when unnamed; global 0); He box lhllc key check requested; class A flips (ck_impl_maxit 24, ck_impl_conserve 1, bbot 0) + sponge_bottom default false -> branch defaults-0930 (agent), merge after gates.

USER 09-30 ~09:45: FRESH 1x and 10x productions WITH the radial fix: w1xf 12037717 (+12037718), w10xf 12037719 (+12037720), inputs w121prod_{1x,10x}_x1f.athinput (lhllc_x1_phi_min = 1), binary athena.gpu.x1f = lhllc-x1-phi 36dd494d; smokes clean (dt 15.5 / 8.4 s). Old w1x/w10x kept as the pre-fix reference.
**DEFAULTS MERGED 09-30 (9bf6db23 + docs 2505ca13, NOTE-2026-09-30-defaults.md):** ck_impl_maxit 24, ck_impl_conserve 1, bbot default 0, sponge_bottom false; gates bitwise with old values named.

09-30 ~09:45: w1xf/w10xf/w3xk binaries swapped (before start) to w121-build-0930 44be5b99 (local build-only merge of rt-integration 2505ca13 + lhllc-x1-phi), md5 ad3e84f7; smoke clean.

**MERGED 09-30 ~10:15: radial fix on rt-integration (aa9d51b7 + NOTE aabe955f), dhj-only default 1, global 0; lhlld floor on chi (bitwise off). Queued w1xf/w10xf/w3xk binary (44be5b99, md5 ad3e84f7) = same src as the merge.
