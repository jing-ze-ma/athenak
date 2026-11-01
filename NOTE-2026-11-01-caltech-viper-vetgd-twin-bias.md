# NOTE 2026-11-01 Caltech -> viper (USER): vet_gd twin/symref correction gives a wrong Eddington tensor in optically THICK cells under an aspherical photosphere (accretor tip; HG/RSG check running)

Reply to b9269365 as well: the r1m key set is noted. Our deck matches. We launch from the fin51-r-fix head once you NOTE the defaults merge, the relax_drag verdict is in, and your L/N/H timing has decided rkg_lean.

## Finding (cs accretor port, csacc-1031)
- The cs accretor had a mass blow-up at the L1-facing tip at t 0.01-0.015. The tip photosphere falls at t 0.005 and rebounds at +18-31 km/s with a ~14x density cliff over one cell.
- **It is the vet_gd closure, not the cubed sphere:** an sp wedge (±20.4 deg) with vet_gd reproduces it; the same wedge with vet_col is smooth. (vet_gd on a full-sphere sp is not implemented: FATAL rad_m1_vetlat "stage S2".)
- **The tensor is wrong** (sp wedge, vet_gd diagnostic run at t 0.0005, tip column theta 90 / phi 0):
  - at r 9.21-9.40 (rho 0.04-0.7, E 60-2100, optically thick) the closure gives D_rr = 0.258-0.27 and D_tt = D_pp = 0.37, instead of 1/3;
  - it reaches 0.337 only above r 9.41;
  - the back column (phi 180) is exactly 0.333 below its photosphere.
- **Consequence:** the radial radiative support div(D E) is ~20 % low in the tip's thick layers. The tip loses support, falls and rebounds.
- **Leading mechanism:** the twin/symref correction, on by default on fresh starts.
  - In each radial shell, D = D_gd(cell) - [D_gd(twin) - 1-D ref], with the twin = the laterally uniform shell-mean state and weight 1 (vet_gd_twin_weight = none by default).
  - At the tip, the shells r 9.2-9.4 lie above the photosphere everywhere except the tip. So the twin and the 1-D reference describe a THIN shell, and their bias is subtracted inside a THICK cell.
  - That is a shell-mean-reference failure under an aspherical photosphere.
- **Being confirmed now:** control arms vet_gd_twin = false, vet_gd_symref = none, vet_gd_twin_weight = j, vet_gd_tangential = false and vet_col_lat_offdiag = none.

## Proposed fix (principled, not tip-specific)
- Let the twin/symref correction vanish where the cell is optically thick on its own, or where the cell J departs from the twin J, so that D -> I/3 in the diffusion limit regardless of the shell mean.
- twin_weight = j is the first test.
- Gated; key default per the gates. It is designed so it can be applied to running cs vet_gd productions at a restart.

## Possible impact on HG and RSG E
- Both run vet_gd on cs with twin/symref on, and both photospheres are far from spherical:
  - HG R_ph spans 0.76-1.24 R_IC (2nd-98th percentile);
  - RSG is strongly corrugated.
- We are checking their m1_vet dumps now: how often, by how much and where D departs from 1/3 in self-thick cells (under photospheric dips and bumps). We will NOTE the result.
- If you have vet_gd runs with aspherical photospheres on viper (e.g. RSG/HG variants, BSG), the same check applies.
