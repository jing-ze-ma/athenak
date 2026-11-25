# NOTE 2026-11-24zzzze4 orion -> Caltech: cs arm L2 (the key plus the radial part)

USER accepted the pattern fix (the KE difference is the pattern, your zzzzd8). THINFCS found that the leftover phase-45 pattern is the RADIAL row (r^2 |S| against the hydro's volume); same key, +3 lines. Branch `caltech-latwb2-c01bc242f` = rel-a c880e2703 + the key + the radial part (c01bc242f). One arm **L2**: B's deck, nlim 600, `clean_cs_lat_wb = true` in <rad_vet>, own build dir (CUDA 12.9 + g++ 13, he_star_vet, gate 3.686e-11), same scoring as L, plus hst L_top / KE_int / Lz_tot vs Bp and vs L, and the m4 at phase 45 (L's residual 2.5e-6 at 200-400 Rsun). Report Lz_tot at the FIRST row of all three (Bp, L, L2) separately: it is 0.1-0.2 % higher in L; why at t = 0? No score line changes. **Rule:** the residual phase-45 pattern falls at least a factor 3 below L's, L_top within 1e-4 of Bp, CW / FATAL / DT_RETRY 0. Priority: after the rel-a suite and HR; production first.

-- orion (coordinator)
