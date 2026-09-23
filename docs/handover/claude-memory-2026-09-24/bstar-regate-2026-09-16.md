---
name: bstar-regate-2026-09-16
description: B-star final_gate re-gate after the bottom-flux fix (real date 09-14/15; the memory index dates are one day ahead) - G3f/G3p/G3q/G3o results, PCR 2.07x on the gate box, tol 1e-6 rejected, PCR solver ALSO lacked the fix (partition.hpp)
metadata:
  type: project
---

bench/bstar_fecz/final_gate, `python3 analyze.py <arms>`; physics columns (v_rms, Mmax, dT_tau, P_nyq) scatter 15 % between arms of this short run - the hard gates are Ftop/Fb, no FATAL, sub/c 49.5, and zone-cycles/cpu_s (the analyze.py ms/cyc of G3 = 115 is from a RESTARTED partial out.txt; use zone-cycles/s: G3 1.04e5 = G3f 1.06e5, so the fix costs nothing).
- G0f, G0p (mode 0): bitwise = G0ref.
- G3f (fix, thomas, wt_rgbox build_hip_botfix): Ftop/Fb 1.00000, clean.
- G3p (fix + PCR nseg32, wt_m3acc build_hip_pcrfix, md5 c9a196b6): Ftop/Fb 1.00000, 2.18e5 zone-cycles/s = 2.07x Thomas on the gate box; BUT the PCR solver path (src/utils/two_stream_column_partition.hpp:115-120) forms dbdtau without the bot_flux override -> on He (Q1p) it collapses bitwise like the unfixed run; G3p passed by the B-star coincidence. Patch: `if (c.bot_flux > 0.0) dbdtau = 3.0*c.bot_flux/(4.0*M_PI);` after line 119. Must be re-gated (G3p2, Q1p2).
- G3q (PCR + rt_impl_tol 1e-6): NO speed gain (98.0 vs 98.6 ms/cyc; Newton already 2-3 passes), P_nyq 1e-3 vs 1e-8 -> REJECTED.
- G3o (PCR + rt_col3_once) REJECTED 09-15 (B star, real date): v_rms 3.0e3 vs 7.5e2 (4x floor), Ftop/Fb 1.0019, only 6 % faster (93.0 vs 98.6 ms/cyc, 2.31e5 vs 2.18e5 zone-cycles/s); the operator-split column brings the pulsation floor back. Earlier note: He Q1o collapsed for the PCR reason, not the once path (once path traced sound: same call, full dt after the last stage via Hydro::RTStrangSplit); B-star G3o ~95 ms/cyc; He Q1o 23.8 vs Q1p 31.3 ms/cycle (1.3x). Needs re-gate after the PCR patch.
- Newton passes cannot be cut further; "Newton = 48 % of the column" came from a TRUNCATING maxit-2 arm.
See [[he-mode3-flux-collapse]], [[bstar-prod-cost-profile]], [[pcr-timing-result]].

TWIN-PATH REVIEW 09-15 (Thomas vs PCR vs sweep, member by member): CLEAN. bvals window latent holes found and fixed in the working tree (build_hip_bvfix): the window must be restricted to same-level neighbours (fine ox1==0 ranges span is..ie on the sender when cnx1==ng while the receiver uses the coarse range) and refused under cubed sphere (FillPanelCornersCC reads the full strips). Both conditions are always true on the B-star box, so bitwise-neutral there. Cosmetic, not done: Dtop[1]/Ucut[1] uninitialised for nq==1 on the Thomas path (unread), stale '12 reduction slots' comment, conduction_transverse passes (is,ie) where (-1,-1) is meant.
MPI GATE: 2 ranks x 2 MeshBlocks, prod9ba8 vs pcrfix2 bitwise (hst + bin payloads). RESTART GATE pending.
G3e (7ca9c57f = restart fix, md5 924480c8): bitwise = G3d until the first rst write (t=1000), then departs by the cache invalidation as designed; v_rms 861, Mmax 0.236, Ftop 1.00000, 2.59e5 zc/s. PASSED on 1 rank. RELEASE BLOCKED only by the decomposition dependence (4 ranks x 1 block).
