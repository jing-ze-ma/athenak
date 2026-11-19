# NOTE 2026-11-19 Orion -> viper (USER): CORRECTION to 14bd22215 (advective upwinding on the staggered lines)

Point 2 of 14bd22215 was wrong. Our F-eliminated staggered line (prod face, no F solve) DOES upwind the advective flux:
on that line rad_m1.cpp reads `cln_upw = cln_split || cln_pform` (commit 1c02017f4, "prod face: the advective donor-cell share
and the implicit gas velocity on"); the `cln_upw = cln_split` I quoted is the older form on the F-solve branch. Its face weights
also take the least donor share that restores al >= 0 >= ar, which is effectively your Peclet > 2 fix.

So BOTH our staggered lines upwind the advective term, and the F-eliminated line's RSG failure (a Picard limit cycle at dt 888.9)
ran WITH that upwinding. Its negative right-hand side comes from the cell row's terms (the work of the implicit gas-velocity kick
and the exchange), not from a centred face advection. Your diagnosis for clean-cc (centred advection at Peclet ~50) therefore does
not explain our F-elim failure, but it is consistent with the F-solve line holding prod's dt (upwinded, and the kick handled inside
the F solve). Points 1, 3 and 4 of 14bd22215 stand.

-- Orion
