# NOTE 2026-10-18 raven: accretor surface + 0.30 test deck, contact eint<=0 burst diagnosed

Follow-up to NOTE-2026-10-17-raven-accretor-sector-light.md (USER-approved 10-17 night).
- Branch accretor-sector-light-1017 @ dfc0a217, PUSHED, NOT merged.
- Full numbers + RESUME: /raven/ptmp/jinma/accbig_1016/NOTE.md "RESULTS 2026-10-18".
- Design: accretor_rhd_design.md section 17 on the branch.

## Code (DEBUG only; gate BITWISE vs 62119e93: stream switch-on + A, hst/bins/cycle lines identical)
- `rad_m1/dbg_eint_cells = N` (ca6f0500): the cells whose implicit gas write-back gives eint <= 0.
- With it, `dbg_t2_admiss` (the E-row-by-term listing) now also runs under time_scheme = be (dfd60ff3).

## (1) TEST deck v4/run/cal30s.athinput = surface + donor_irr_scale 0.30 (old decks untouched)
- Static and stream to t 0.0268: dt 7.02e-6 flat, 0 retries.
- To t 0.036 only with the FALLBACK dt_min 1e-8:
  - 49 eint<=0 retries at contact (t 0.0270-0.0273), then dt 1.5-1.65e-6 (v3 0.84-1.6e-6);
  - NC 0, Mnzf 19.7958;
  - thin-atmosphere max|v| 93-224 km/s before contact, 465-527 after (v3 376-1117);
  - impact = v3 (heating from t 0.032, phi 13.8 -> 15.0 deg, e/rho 4039 at t 0.036).

## (2) Contact burst: cause
**Failing cells:** the stream's contact front, phi 6-10 deg, r 9.4-9.85, T 3.6-6 kK. In 90 % of them the
SOLVED E is negative.
1. The negative E is floored.
2. The exchange becomes pure emission, which is about the cell's eint.
3. The write-back gives eint < 0 and dt_retry fires.

**Why E < 0:** the E row is a small residual of two terms ~10^3 x E. The old face fluxes theta F^n
(transverse -660..+330, radial +650..-320) appear with theta ~ 1 in thin cells, at c dt/dx ~ 190 (r) / 16 (phi),
and the transverse faces have no free-streaming bound. When the cold opaque front moves into the floor gap the
stored fluxes no longer balance.

**Why only with dim sector light:** the stream arrives cold (e/rho 1.55e3 vs 2.3e3 in v3/scaled), in the
recombination range of the opacity table, into a field 5-20x dimmer.

**Ruled out:**
- `implicit_g0_exchange = false`: FATAL.
- `implicit_gas_newton_u4 = true`: FATAL.

## Proposed fix (USER: change to the shared M1 solver, opt-in, IDORT validation)
A positivity-preserving old-face-flux term, by either:
- (a) `implicit_trans_limit = lp` (LP face limiter) extended to the sp wedge, plus radial faces;
- (b) FCT-type scaling of the outgoing theta F^n.

Interim fallback: dt_min 1e-8.

## (3) Cleanup
~38 GB of restart files deleted (list in NOTE.md). Kept:
- surface/rst/00004
- surface_imp8/rst/00005
- cal30s/static/rst/00001
- cal30s/stream/rst/00005 (contact-study start)
- cal30s/stream_fb/rst/00010
