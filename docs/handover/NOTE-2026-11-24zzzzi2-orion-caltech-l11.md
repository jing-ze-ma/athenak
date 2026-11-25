# NOTE 2026-11-24zzzzi2 orion -> Caltech: bsg leak: add arm L11 (prime suspect), run it FIRST after LZB0f

Suspect (LZLEAK, code reading, not measured): the take-back of the x1 wall / top-inflow fluxes (star_zflux, problem/star_wall_zero_flux = true) removes the momentum with the cell's lever arm, while the flux was given with the face's arm: leftover about dr/(2r) ~ 4e-3 of the taken-back L_z flux, biased, independent of floors / FOFC / merge / tolerance; r / |Lz_hflx dt| was 4.6e-3 mean.
| arm | change |
|---|---|
| L11 | problem/star_wall_zero_flux=false (no take-back; +100 from bsg3d.00004.rst; r and Lz_hflx as before). If the code FATALs (a wall needs it), say so and skip |

Expected if the suspect is right: r falls by orders of magnitude, although the wall fluxes are then real fluxes (mass and energy through the walls change too: diagnostic only, not a proposal). If r stays 3e33, the suspect is out. Also report M_tot and Etot drift vs the booked fluxes for L11 so we know that the arm is consistent.
The rest of the battery (L1-L10, zzzzi1) stays, after L11.

-- orion (coordinator)
