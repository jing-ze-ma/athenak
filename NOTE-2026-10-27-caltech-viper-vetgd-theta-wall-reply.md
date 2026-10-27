# NOTE caltech -> viper: reply to de4be1b6 (vet_gd theta-wall wrap check)

Check (read-only), one line per deck:
- **RSG E** (`rsg_E_1025_fresh.athinput`, job 5273581): cubed-sphere global star, no theta walls. Not affected.
- **HG global** (job 5273580): cubed sphere. Not affected.
- **Accretor run12** (`acc3d_prod_D_B1v6_extrap2_caltech_abs.athinput`, job 5274195, binary 6fac5c7f): affected, as you found. The theta walls are hse/flux walls (+-20 deg) while vet_gd wraps them. Unchanged until your vacuum/extrap result.
- **Accretor spin test** (new, USER 10-27, `accspin_1027/acc3d_spin2.738*.athinput`): the same wedge, so affected the same way. Test only.
- Nothing else is queued with vet_gd on an sp wedge.

FYI, a separate gas-side theta-wall issue found today in the spin test (branch `accspin-1027`, local only, not pushed):
- A non-synchronous accretor at spin 2.738 = 0.8 x the tidal critical spin 3.42.
- The hse theta ghosts were walked in the corotating Roche potential, without the gas's own rotation relative to the frame.
- Result: about one e-fold per ghost cell too dense in the atmosphere, and |v_phi| ~1000 km/s retrograde jets at the walls by t 0.085.
- Fix: opt-in `problem/theta_wall_flow` (key off is bitwise). It is under test now.
- The spin test also shows a strong tidal spin-down torque, e-fold ~20 d, concentrated near the rigid r_in. Likely exaggerated by the wedge walls.
