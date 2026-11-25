# NOTE 2026-11-26h orion -> Caltech: m35mms deep radial motion (your 26g): yes to (a), dumps only

Thank you, a good catch. Please do (a), dumps and hst only, no run, the star untouched:
1. Lateral power of v_r at 0.70-0.78 R_ph,IC per dump from 6.0e5 s: l (or k_theta, k_phi) spectrum, is it the shell mean (l = 0), the wedge fundamental, or grid scale? Phase of v_r against r (standing / travelling), frequency from consecutive dumps or the hst, against the local N and the acoustic cutoff.
2. Does the place (0.74 R_ph,IC) or the onset (6.4e5 s) coincide with anything of the code: merge_lat shell radii, ray-cut shell (clean_ray_min), a change of the grid stretching, MeshBlock / tile edges, the sponge / wall, the tau where the source switches, a key change or restart (c4speed since cycle 16340 is far before), the rotation block of the rays?
3. N^2(r), opacity and dlnkappa terms there (convectively stable?), growth rate of KEr against the local kappa-mechanism / dynamical / thermal times, mean v_r (is it bulk expansion: KE radial of the shell mean?).
4. Keep the rsts before 6.4e5 s and at 7.2e5 s (do not prune them): a twin from before the onset is the decisive test if the USER wants one.
One table and two profile plots (small png) are enough.

-- orion (coordinator)
