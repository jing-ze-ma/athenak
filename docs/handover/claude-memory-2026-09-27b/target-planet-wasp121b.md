---
name: target-planet-wasp121b
description: User decision 09-25 -- the dhj model targets WASP-121b (all planet/star/rotation parameters from literature); R_p placed at the transit (chord tau ~0.56) pressure computed with our model
metadata:
  type: project
---
User 09-25: "match all parameters to wasp121b; decide based on literature where to put the observed Rp".
Old model planet: Teq 2500 K, 0.66 M_J, 1.42 R_J at 1 bar, g 942 at r_in (not a real planet). WASP-121b (quick
search, to be verified by the setup agent): M 1.157 M_J, R 1.753 R_J, Teq ~2358 K (Delrez+2016), P_orb 1.27 d
(tidally locked -> P_rot), g ~930 cm/s^2. Setup agent: /viper/ptmp2/jinma/wasp121_0925, branch planet-setup
(planet_setup.py: M_p, R_p at p_ref -> x1min/x1max/ap/grav + redesigned grids; new WASP-121b RCE ic_profile).

**How to apply:** all new dhj inputs use the WASP-121b parameters from that README; the sponge tests should also run
on WASP-121b (ask if unsure). Related: [[sparc-sponge-campaign-0925]], [[dhj-deep-convective-verdict]].

**Albedo (user 09-25): use the OBSERVED Bond albedo for WASP-121b** (Frazier+2026 A_B 0.277 +- 0.016 candidate; setup agent picks with reasons) via new problem/albedo overriding the Parmentier+2015 get_albedo fit (~0.02); RCE ic_profile recomputed with the same (1-A); Teq stays literature zero-albedo.

**WASP-121b setup DONE 09-25 (/viper/ptmp2/jinma/wasp121_0925/README.md; merged dhj-albedo + planet-setup 038148c6):**
Sing+2024 (NASA default, homogeneous): M 1.170 M_J, R 1.742 R_J (JWST NIRSpec NRS1), a 0.02571 AU, Teff 6628 K,
R* 1.461 Rsun, P 1.27492504 d (tidally locked), Teq 2409.3 K. Atmosphere 10x solar (met 1.0, eos_xh/yhe 0.6588/0.2218,
upstream Exo-FMS 10x ck + CE + generated 10x hiT). Stellar SED: WASP-121 band-flux file (ck_star_teff = 0). Bond
albedo 0.277 (Splinter+2025 JWST NIRISS phase curve). p_ref = 1.8e-4 bar (NRS1 transit radius from our ck on the
10x RCE; 1x would be 7e-4). Keys: x1min = ap 1.116018e10, grav 1190.073, omega 5.704026e-5, x1max 1.530921e10,
dfloor 1e-16. Grids: SPARC nx1 74 (dt ~14 s), production nx1 256 (~10 cells/H, dt ~4.6 s). run.sub needs INP,
grid_w121.env ($G), PR 1.101535e5, output dts 1.10154e3 / 5.50768e4. Open: IC table only to 283 bar vs equatorial wall
~630 bar (bulge) -> being extended; binary must be rebuilt with dhj-albedo.
