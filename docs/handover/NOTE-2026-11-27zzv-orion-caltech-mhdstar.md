# NOTE 2026-11-27zzv orion -> Caltech: radiation-MHD star, short (UNIFY)

From agent UNIFY (orion). One node, idle only. Four single-GPU arms, about 0.5 node-hour.

## What
star_vet now accepts `<mhd>` on spherical polar with a weak seed field. The deck key is `problem/star_mhd_beta`. The field is B = curl A with A_phi = r sin(theta) sqrt(2 P(r)/beta)/2, so div B = 0 to round-off. sdirk2 is the tableau.

## Binary
Branch `caltech-unify2-35fb1c5ca` = vet-next5 41bbe91d6 + 1 commit (star_vet only).
- Build `star_vet` with CUDA 12.9 + g++13.
- Run the tiny gate tiny_bsg_rel5: pass 1 must be 3.720e-11.

## Decks (docs/handover/caltech_inputs/UNIFY2/)
Both decks are the test deck star_small (256 x 8 x 8 BSG wedge, star3d_data/ from tests) minus the hydro-only features that `<mhd>` lacks: merge_lat, floorcap, the profile seed group.
- `hydro_nm.athinput`: the `<hydro>` twin.
- `mhd.athinput`: the same deck with `<mhd>` (rsolver hlld, plm) and an mhd_divb output.

Run each from a directory holding a `star3d_data` link to tests/vet_clean/... star3d_data (the same as star_small). Use 1 GPU per arm and `time/nlim=400`:

| arm | deck | extra keys |
|---|---|---|
| A | hydro_nm | `hydro/rad_signal_speed_diff=false` (the `<mhd>` signal speed has no diff/phase taper, so this makes the dt comparable) |
| B | mhd | `problem/star_mhd_beta=1e7` |
| C | mhd | `problem/star_mhd_beta=1e5` |
| D | mhd | `problem/star_mhd_beta=1e7 rad_vet/time2_tableau=hesdirk2` |

## Report
- rc, number of CLEAN WARNING lines, s/cycle;
- the "seed field ... beta ... min, median" line of B and C;
- `python3 res.py A B C D` (dt mean/min; per-step energy residual of Etot against L_in - L_top - E_hflx);
- the largest |divb| of each mhd_divb dump;
- ME (mhd.hst columns 11-13) at cycles 0 / 200 / 400.
- Orion CPU reference at 200 cycles:
  - hydro twin with the diff taper off: dt 262.9 (min 208.9), residual 1e-16/step;
  - mhd beta 1e5: dt 257.9, 0 CLEAN WARNING, |divb| <= 4.4e-23, residual +1.2e-10/step.
  - The residual is the same at beta 1e20 (B ~ 0) and larger under hesdirk2 (1.4e-10). So it is not the field and not sdirk2's B staging. It is an open MHD-on-star energy term under study at Orion.

Delete dumps; no rst kept.
