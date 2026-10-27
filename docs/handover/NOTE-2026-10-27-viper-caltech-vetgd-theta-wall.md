# NOTE viper -> caltech (USER): vet_gd formal solution is theta-PERIODIC at every theta wall; this affects the accretor (run12 family)

FYI. Nothing to change in the running jobs yet. Please reply with your answer to the check at the end.

## The bug
- In the sp wedge, the `vet_gd` formal solution treats every theta wall as **periodic**, whatever the M1/gas boundary is.
- Code at 6fac5c7ff, `src/rad_m1/rad_m1_vetgd.cpp`:
  - line 1249, `VgdTabCompute`: `if (th > t1) {thw = th - (t1 - t0);}`;
  - line 1677, `VetGdHaloInit`: `l2 = (l.lx2 + dj + nb2) % nb2`.
- So rays near a theta wall read the radiation of the OPPOSITE wall.

## Where it was found
- A thin clump on a reflecting equator wall, on the CPU, compared with the mirrored domain (no wall), which is exact:
  - formal-solution J in the wall layer: 14-16 % low;
  - moment E in the wall layer: +86 / +85 / +152 % (median) at nx1 32 / 64 / 128, so it does not converge;
  - E-floor clips: 0 / 1004 / 121309.
- The moment solve, the M1 ghost fill and the geometric terms are clean. The formal-solution wrap is the whole effect.

## Who is affected
- **Accretor wedge family** (run12 `acc3d_prod_D_B1v6_extrap2_*`, and the hold/vphi/pfb/csfix accretor decks). The theta walls at ±20 deg are FLUX walls for gas and M1, while the formal solution wraps them.
- **Size** (viper, X6 restart, 30 cycles, production keys): switching the wall treatment changes E in the wall rows by median 0.5 %, p90 10-11 %, max 41-44 %. In the outer thin region of the wall rows the median is 4.8 %; the mid row changes by 0.015 %; gas density by at most 2 %. This was still growing at 30 cycles, so the steady-state size is unknown.
- **Not affected:** the HG wedge (theta-periodic by design, so the wrap is correct) and HG global (cubed sphere).
- RSG E was not checked on viper; see the check at the end.

## Fix status (fork only, opt-in; default = old behaviour, bitwise)
- `fork/wallbug-1027` = b104357b5 adds `<rad_m1>/vet_gd_theta_wall = periodic | mirror`. Key off is bitwise on the accretor deck, at +0.7 % cost (noise).
- **Do NOT use `mirror` for the accretor.** It is correct for a REFLECTING wall only, and the accretor walls are outflow walls. With the production GL10x12 staggered direction set it is also only approximate (p90 error 0.10 -> 0.43 with resolution); it is exact only with an odd `vet_gd_gl_nmu` or `vet_gd_gl_stagger = false`.
- **In progress on viper:** `vet_gd_theta_wall = vacuum`, which gives zero incoming intensity at the wall, plus possibly `extrap`. It will be measured against periodic on the run12 restart. A follow-up NOTE or TASK will say whether run12 should switch at a restart. No binary or deck change until then.

## Check, please (read-only)
Do any of your running or queued production decks combine `vet_gd = true` on an sp wedge whose theta boundaries (`ix2_bc` / `ox2_bc`, or the `<rad_m1>` implicit_bc_x2*) are NOT periodic, other than the accretor? RSG E in particular. A one-line answer per deck is enough.
