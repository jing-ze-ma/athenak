# NOTE 2026-11-28a orion -> Caltech: radiation-MHD star A/B, ported features (UNIFY)

From agent UNIFY (orion). Short run: idle node only, 4 single-GPU arms (sp) plus one cs pair if a node is free, about 0.5 node-hour.

## Binary
Branch `caltech-unify3-eef7f8177`, which is e02e092e8 + 1 commit. That commit adds, under `<mhd>`:
- dfloor_flux_limit, floorcap, dt_floor_courant and the signal-speed taper;
- the star_vet <mhd> features on sp and cs;
- Strang RKG with <rad_vet>.

Build `star_vet` with CUDA 12.9 + g++13. The tiny gate is **tiny_bsg_rel6** (time_scheme = sdirk2); pass 1 must be 3.720e-11.

## Decks (docs/handover/caltech_inputs/UNIFY3/; run with a star3d_data link as before)
| arm | deck | ranks | keys |
|---|---|---|---|
| A | hydro2 (star_small without merge_lat and the seed group) | 1 | time/nlim=400 |
| B | mhd2 (the same deck with `<mhd>`, hlld, star_mhd_beta 1e5 = the smallest beta) | 1 | time/nlim=400 |
| C | mhd2 | 1 | time/nlim=400 problem/star_mhd_beta=1e7 |
| D | mhd2 | 1 | time/nlim=400 mhd/rsolver=hlle |
| E/F (optional) | hcs / mcs (star3d_cs without merge_lat) | 6 | `mesh/nx1=274 meshblock/nx1=274 mesh/f_stretch_r_head_nin=296 mesh/f_stretch_r_head_ij=55 mesh/nx2=16 mesh/nx3=16 meshblock/nx2=16 meshblock/nx3=16 rad_vet/vet_gd_allow_clamp=true time/nlim=40` |

## Report
- rc, CLEAN WARNING count, s/cycle;
- `python3 res.py A B C D` (residual per step against L_in - L_top - E_hflx + P_esrc);
- the max |divb| of the mhd_divb dumps;
- for E/F, the residuals from the *.user.hst files.

Orion CPU reference at 200 cycles:
- A: 9e-18 per step.
- B: +1.7e-14 per step. It scales with B^2: C gives +1.7e-16. With hlle (D) it is -2.6e-15.
- B and A have the same dt (314.5).
- E and F have the same residual as each other (3.3e-12; the budget misses a cs term).

Delete dumps; no rst kept.
