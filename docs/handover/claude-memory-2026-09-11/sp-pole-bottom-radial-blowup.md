---
name: sp-pole-bottom-radial-blowup
description: sp_mhd_prod (18f5dd21, the first sp MHD production run with the 09-05/06 pole+resistivity+wall fixes) grows a RADIAL field in the two polar rows at the BOTTOM (i<=16), odd azimuthal m 5-9, e-folding every ~0.4 rot from rot 2, 4.9 kG by rot 6, NaN at rot 10.3; sp_mhd_diss (3147de4f) never did. BISECTED to 2a64e2c7..18f5dd21 (2a64e2c7 CLEAN at rot 5, ideal+noavg on 18f5dd21 DIRTY, resistivity EXONERATED); round 2 arms wall_r (4990eb41 = without fefe6a17) and x3shift_r (863e8337) in flight 09-07 05:32
metadata:
  type: project
---

2026-09-06 evening. **sp_mhd_prod (job 11501024, binary 18f5dd21, input = sp_mhd_diss's)
went NaN at rot 10.3** and ran on with NaNs to rot 41.5 (the scancel was blocked for me;
the user must cancel 11501024/11501025). Restarts: rst 00000 = t 0, 00001+ = rot 10+
(NaN) -- no usable early restart. sp_mhd_diss (3147de4f, same input) is clean to 50+ rot.

**Signature (bins every 2 rot, scripts $S/sploc.py, $S/polerow.py, $S/inflow_gate.py).**
History ME1 (radial magnetic energy): diss 8e31 flat; prod 8e31 at rot 3.7 -> 9.7e33 at
4.7 -> 8.4e34 at 5.7 -> 1.0e35 saturating. Polar-row max|B_r| at i=0, j=0 (north; south
mirror-identical): rot 2: 5 G (diss 11 at rot 6), rot 4: 171 G, rot 6: 4840 G; j=1 3668,
j=2 1702, j=3 838, j=5 10 -> confined to the 4 polar-most rows; by radius max|B| 4.9 kG at
i=0 falling to 2.4 kG at i=7, ~100 G at i=16, nothing at i=48. The phi-spectrum of B_r in
the pole row is ODD m = 7, 5, 9 (and 21 later), NOT the Nyquist checkerboard of
[[sp-polar-field-blowup]] and not the m=1 dipole that diss shows. NO bottom inflow (min
v_r -0.09 km/s at i=0 through rot 10) -- unlike the cs defect [[cs-mhd-bottom-inflow]].
B_theta, B_phi in the same cells are 5-10x smaller than B_r. Floors/C2P counters silent.

**Commits between the binaries (3147de4f..18f5dd21) that touch sp MHD at default
settings:** 2a64e2c7 (3D curl: J_r, J_phi regain their theta derivatives -- see
[[resistivity-3d-curl-missing-terms]]), 3694cc40 (pole dual edge area), de667f32
(polar_emf_diss third differences, default ON), b84a0502 (use_polar_average_eresist
default ON), a8cfb83e + 863e8337 (x3-face states in the local basis / at the face
midpoint, ALL variables, no flag), b61c5d67 + 581a1e43 (pole-edge J_r sign / cap sector),
fefe6a17 (reflecting-wall mirror). Flags that exist: mhd/polar_emf_diss,
mesh/use_polar_average_eresist, mesh/polar_quadratic_recon, mhd/sp_cart_{polar,all}_momentum.

**In flight (submitted ~20:00-20:50, 3.5 h apu1 each, 5 rot from scratch, bins 1/rot):**
bench/sp_pole_bisect/{nodiss 11516338, noavg 11517118 (11516339 FAILED at parse: use_polar_average_eresist is not in the input, so the cmdline override is refused -- it is now IN the input file), ideal 11516341, cfl 11516450 (cfl 0.15)}
on the 18f5dd21 binary, and curl2a64 11516469 on the 2a64e2c7 binary (BUILT, in
bench/wt_2a64e2c7/build_gpu/src/athena; worktree has kokkos symlinked to the main repo). Gate: history ME1 at rot 4.7 (prod 9.7e33 vs diss 8.3e31) and polar-row |B_r| at
i=0 in the rot-4 dump (prod 171 G, diss 19 G): `python3 $S/polerow.py <bin>`. ideal clean
=> resistive path (then run the 2a64e2c7 binary being built in
bench/wt_2a64e2c7/build_gpu, log $S/build_2a64.log, kokkos symlinked); ideal dirty => CT /
x3-face-shift / wall (863e8337, a8cfb83e, fefe6a17; no flags -> binary bisection).
Hypotheses, unranked: the new d(B_r)/dtheta term of J_phi across the polar ghost row
(ghost x1f exchange sign?), the cap-sector pole-edge J_r, resistive dt at the stretched
polar row (8.1 deg cells). Why the BOTTOM: unknown (deep = hot = low eta, so not eta).

**00:10 09-07: the 3.5 h arms were TOO SHORT** -- ~1 rot/h early (dt 2-3 s), none reached
rot 4.7 (ideal 3.13, nodiss 2.29 but ME1 1.4e34 already = blows up FASTER without the polar
EMF dissipation, so that is a stabiliser not the cause; noavg 2.0, curl2a64 2.5, cfl 1.2 all
~6-7e31 = normal). No restarts were written. RESUBMITTED as {ideal,noavg,cfl,curl2a64}_r
(jobs 11520726-9, 7:45 h, rst every rot, submit.sh resumes from the newest rst) -- if
maintenance (09-07 12:00) cuts them, resubmit the same script after 09-12 12:00.


## 05:30 09-07: ROUND 1 RESULT. 2a64e2c7 is CLEAN, 18f5dd21 is DIRTY, resistivity is innocent

History ME1 at rot 4.7 / 5.0: ideal_r 1.6e34 / 4.7e34, noavg_r 8.5e33 / 3.9e34,
curl2a64_r 8.8e31 / 7.9e31 (flat, = diss). Polar-row |B_r| at i=0 in the rot-4 dump
(polerow.py): ideal_r rms 48 G top m [7,5,9], noavg_r 54 G, curl2a64_r 5.9 G with
m=1 (the tilt). cfl_r (cfl 0.15) only reached rot 1.9 (ME1 7.4e31, normal; dt 0.5-0.8 s,
too slow to gate; cut by the maintenance, resumable from rst). ideal_r has NO
resistivity and blows up => the culprit is in the ideal-MHD sp commits of
2a64e2c7..18f5dd21 (chronological): 3694cc40 (resistive), a8cfb83e (x3-face states in the
local basis), de667f32 (polar_emf_diss 3rd diff), b84a0502 (resistive), eb8e3cb1
(Cartesian polar-row update triads), 863e8337 (x3-face shift, all variables), 0f9959ee
(sp_cart_all_momentum, flag), b61c5d67/581a1e43 (resistive), fefe6a17 (REFLECTING WALL
mirror -- the blow-up is at the BOTTOM wall), 87c18855/18f5dd21 (tests/style).

ROUND 2 (submitted 05:32, 6 h, sized to end before the 12:00 maintenance; ideal_r's input,
5 rot, rst every rot, submit.sh resumes): wall_r = 4990eb41 = 18f5dd21 minus fefe6a17
(job 11526902, bench/sp_pole_bisect/wall_r, worktree bench/wt_4990eb41);
x3shift_r = 863e8337 (job 11526903, bench/wt_863e8337). Read: wall_r clean => fefe6a17
(the wall mirror) is the cause. wall_r dirty & x3shift_r clean => 0f9959ee. Both dirty
=> one of a8cfb83e/de667f32/eb8e3cb1/863e8337 (one more round). Gate as before: ME1 at
rot 4.7 and polerow.py on bin 00004. Trap: waiting on a build log with grep "rror"
matches cmake's own output -- wait on BUILD_DONE only.
