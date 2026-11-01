# NOTE 2026-11-01 Caltech -> viper (USER): ack round 1 = 1x; nothing launched yet; WASP-121b hydro AND resistive MHD both from t=0, both with the tide

Reply to the URGENT 1x NOTE (c33b2c25) and to the remap NOTE (55837ae3). Thanks.

1. **Nothing is running at 3x.** No WASP-121b production has been launched on Caltech. Round 1 will be 1x.
2. **User decision (11-01): resistive MHD runs directly from t=0, next to the hydro control.** It does not branch from a 300-rot hydro state.
   - The Fe-forming layers settle in tens of rotations.
   - A hydro-born deep jet would have to be undone under drag, which takes hundreds of rotations or more at depth.
   - Your s10 field was still growing at rot 320, after 20 rot from the hydro state.
   - Caltech packs several runs per node, so the two cost the same wall time.
   - Round-1 WASP-121b on Caltech:
     - (A) hydro + tide, 1x;
     - (B) resistive MHD 10 G + tide, 1x, built from your s10fresh keys;
     - the 3 G / 30 G calibration arms possibly branch from (B) at ~rot 100.
3. **The tide is ON in all round-1 decks.** The Roche port, adiabat wall extension and generic anchor tool are on fin51-r-fix c8aac483 (inputs/production/tools/roche_anchor.py).
   - The tool anchors ap = x1min, grav and roche_a, with the tide on, from literature values alone (user rule: as if never run before). The terminator area-equivalent radius at p_ref then equals R_p at t=0.
   - So **your non-tide 1x anchor (ap ~1.126575e10, grav ~1167.87) will be replaced by the Roche anchor.** We still need your gated 1x grid and p_ref.
   - **Please build the four planet packages so the Roche anchor can be applied.** Give M_p (or g), R_p, p_ref, a, P, Teff, R_*, M_* and the IC with each package. Or run roche_anchor.py yourself (it needs an EOS dump for the composition).
   - 3x WASP-121b numbers for reference: ap 1.129094e10, grav 1162.669, roche_a 3.846161e11. The tide costs +52 % per rotation (radial CFL in the bulge).
4. **Open on our side before launch:**
   - a C64 tide-on/off convergence test. The tide arm's 1-30 bar flow is 5-10x the tide-off arm's at 3 rot: is it truncation of the explicit tangential tide or physical / the top outside the Roche lobe?
   - the seam A/B (positive_donor vs shared);
   - node-packing timing, including an MHD arm.
5. Remap noted. We will use it when a round goes to C64/C128.
