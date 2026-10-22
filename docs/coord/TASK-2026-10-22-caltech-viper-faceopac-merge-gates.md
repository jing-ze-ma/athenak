# TASK 2026-10-22 caltech -> viper (USER): gate + fast-forward fin51-r-fix to faceopac-prod-1022 yourselves (re ea2a7ea4)

USER 10-22: Caltech runs no tests inside its production allocations; viper runs the merge gates. Caltech's gate attempts were
stopped before a verdict (no results to use).

Branch (pushed): fork/faceopac-prod-1022 = b8f75bfa
  c92d5b03 (fin51-r-fix) + de0ee02b kirch-1022 (implicit_face_opac = arith default | kirchhoff)
                         + d33f5f59 harm-1022 (= harmonic)
                         + b8f75bfa: the gsdrag-1021 cs2_dacap drag sites (rad_m1_impl_wb_kernel.hpp ~509/515, added on fin51 after
                           kirch-1022 branched) also use M1KtX1 (same face rho kappa_T as the deposit); key-off identical call.
  One cherry-pick conflict (rad_m1.hpp declarations next to implicit_vimp_kill): both blocks kept.
  NOT included: mlt_fsub_smooth, opac_lag_passes = all, implicit_opac_relax (diagnostic dead ends).
Caltech H200 builds (for reference): c92d5b03 md5 b743d0bc, b8f75bfa md5 0ccfb257.

## Asked of viper
1. Gates on your machine: key-off b8f75bfa vs c92d5b03 bitwise (he + accretor decks, same environment); key-on harmonic C2 == C and
   R == C bitwise (HG or RSG E deck). Your earlier d3df1401 key-on gates PASSED; b8f75bfa adds only the cs2_dacap site change.
2. If PASS: fast-forward fin51-r-fix to b8f75bfa, then rebase + gate prod-extrap2-1022 on top as you proposed; NOTE the hashes.
3. Accretor / BSG harmonic A/B (TASK 9c432f17) continues; default flip only after it.
Also: USER decision implicit_cs2_dacap stays OFF (845f5e49).
