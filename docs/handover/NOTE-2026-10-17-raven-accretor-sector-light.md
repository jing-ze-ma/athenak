# NOTE 2026-10-17 raven: accretor back-sector light switch + confirming tests + vet_col lateral-beam diagnosis

Follow-up to NOTE-2026-10-17-raven-accretor-v4.md (USER-approved 10-17).
- Branch accretor-sector-light-1017 @ 5eeb137e (code 62119e93), off accretor-bigdom-1015 99adf52d. PUSHED, NOT merged.
- Full numbers + RESUME: /raven/ptmp/jinma/accbig_1016/NOTE.md "RESULTS 2026-10-17b".
- Design: docs/dev/accretor_rhd_design.md section 16 on the branch.
- Everything ran on Raven A100 (one platform).

## Code
`problem/stream_sector_light` sets the radiation held in the back-sink sector. It is read only when named.
- `full` (default = v4)
- `scaled`: x donor_irr_scale. This is an upper bound for the stream; its own 3-D factor is ~0.07.
- `surface`: the lobe's donor-light rule, i.e. the v3 rule.
- `dark`: no light held.

**Gate BITWISE** vs v4 5b0d4bc1, key absent, 20 cycles:
- the stream switch-on from the v4 static rst, which exercises the sector: hst, bins and cycle/dt identical;
- v4 gates A and G2: identical.

## Confirming tests (from v4/cal/static rst, s = 0.22)
| arm | dt to contact (t 0.0095-0.0265) | contact t 0.027 | dt after | thin-atm max\|v\| (t 0.013/0.020/0.031) | Ltop (t 0.0107/0.0285) |
|---|---|---|---|---|---|
| v3 | 5.6e-6 -> 1.5e-6 | smooth | 1.2-1.6e-6 | 376/717/1117 | 3.7e5/7.8e5 |
| v4 full | 2.9e-6 -> 3.1e-7 | 49 retries over the run | 4.5-5.5e-7 | 660/1429/2394 | 1.2e6/4.1e6 |
| scaled | 6.7e-6 -> 2.9e-6 | 1 retry | 1.5-1.6e-6 | 193/468/1068 | 3.0e5/7.5e5 |
| surface | 7.02e-6 flat | eint<=0 burst: FATAL at dt_min 1e-7; recovered at 1e-8 (61 retries) | 1.56-1.65e-6 | 80/139/469 | 6.9e4/1.0e5 |
| dark | 7.02e-6 flat | same (FATAL at 1e-7; recovered at 1e-8) | 1.57-1.63e-6 | 77/128/474 | 6.3e4/8.9e4 |

- The impact matches v3 in every arm: contact t 0.0272-0.0284; heating from t ~0.031 at phi 13.6-13.8 (v3 14.0).
- Floor-gas e/rho stays at ~2e3 (v4 full: 1.8e6).
- **The unscaled sector light was the v4 dt killer (confirmed).**
- **New: a contact transient in surface/dark.** It needs dt_min < 1e-7 or a fix; the cause is not located.

**s = 0.30 static round:** facing donor flux at r 10.0 / 10.4 = 1.24x / 1.00x the 3-D target (s = 0.22 gave
0.73 / 0.52).

## Recommendation (USER decides)
**`stream_sector_light = surface` + `donor_irr_scale = 0.30`.**
- It is one rule with the lobe.
- It is physical: the stream's core light is trapped in the thick stream; it must not stream out of thin floor cells.
- It needs a policy for the contact transient: dt_min 1e-8, or locate the cell.
- Fallback: `scaled`. It is v3-like, but its pre-impact dt dips to 2.9e-6, the atmosphere reaches ~1000 km/s, and
  the light is unphysical.

## vet_col lateral beam: diagnosis
1. The tensor is column-local and laterally homogeneous: diag(f_K, (1-f_K)/2, (1-f_K)/2) about r_hat. Lateral
   donor light enters only through E in S, as if it were isotropic local emission.
2. D_r,phi = 0. The sp face equations carry only the diagonal curvature (M1SphCurv, rad_m1_implicit.hpp:906), so
   nothing turns F_phi into F_r. F_phi only feels -(c^2/r) d_phi((1-f_K)/2 E) - c chi F_phi: a lateral
   wave/diffusion along each shell that thin gas barely damps. Light passing the limb circulates round the star
   instead of leaving along straight lines, so there is no shadow.
3. D does not depend on F: the donor-part |dF_phi|/(c dE) reaches 1.0-1.93.

**Prototype, the existing opt-in `vet_col_lat = true`** (supplies D_r,lat): UNUSABLE on the 4-cell theta band.
E falls to e_floor at r 9.6-10 by t 0.0012, with retries from cycle 76.

**Fix options (method change, USER):**
- (A) ray-traced external donor irradiation with exp(-tau) and an exact shadow, the lobe light not M1-transported;
- (B) an F-aligned thin-region closure blend (the thinclos route; it needs debugging first);
- (C) a 3-D wedge nx2 >~ 64 for vet_gd.

No vet_col default changed.
