---
name: red-giant-seed-across-rcb-proven
description: "PROVEN 2026-09-09 (R2a/R2b/R2c vs R0f): the RCB superadiabatic dipole -- and hence the deep 10x-MLT startup transient -- is caused by the vpert SEED crossing the radiative-convective boundary. problem/vpert_rmin = 1.6e12 removes it entirely (|sup| < 3e-5 flat) with or without MLT; quadrature alone changes nothing. Also: ck_nquad=2 killed BOTH no-MLT runs at t=1.44e6 (bisect before production); MLT v3's outer gap is the chi limiter."
metadata:
  type: project
---

| run | seed | nquad | MLT | sup(6) at 1.4e6 | L_out/L | T(tau=2/3) | fate |
|---|---|---|---|---|---|---|---|
| R0f | everywhere | 1 | off | +4.35e-4 (linear) | 0.36-0.60 | 4279 K | fine to 2.1e6 |
| R2c | everywhere | 2 | off | +4.35e-4 (identical) | 0.37-0.68 | 4344 K | **NaN at 1.4417e6** |
| R2a | r > 0.5 R | 2 | off | **-2.3e-6 (flat)** | 0.37-0.66 | 4130 K | **dt collapse 1.4472e6** |
| R2b | r > 0.5 R | 2 | v3 + F_2s | **-4.0e-6 (flat)** | 0.7-1.3 | 6956 K | alive, min dt 17.9 |

- The seed, `v1 = vpert*cs*sin(3*2pi*x/block)`, was imposed across the RCB entropy step at
  r ~ 0.04 R; a standing radial flow through an entropy step advects entropy at a constant
  rate -> linear-in-time dipole from t = 0, blind to any closure.  Cutting the seed below
  0.5 R (5-cell raised-cosine ramp) removes it.  This supersedes the "pile-up" mechanism in
  [[red-giant-deep-onset-is-rcb-pileup]] as the PRIMARY cause: the flux mismatch that
  closure v3 fixed to 1e-4 was never what fed the dipole.
- **ck_nquad = 2 is suspect**: both no-MLT runs with it died within 6000 s of each other
  at t = 1.44e6 (R2c: 46856 non-finite cells; R2a: hydro dt 3e-14), where R0f (nquad 1)
  sailed through.  Being bisected from the last dumps before adopting it in production
  despite its correct diffusion limit ([[red-giant-quadrature-seed-binary]]).
- **MLT v3 gap**: at 1.4e6 R2b carries F_mlt = 1.00 F_req from 0.15 to 0.6 R, but from
  i ~ 160 to 300 F_mlt falls 0.65 -> 0.16 with nothing else carrying (F_2s is exactly 0
  below its icut ~ 285), so the layers above heat and radiate 1.3 L.  Cause: the chi
  limiter evaluated at x_eff = x_thr = 1e-4 with f = D -> chi enormous -> throttled.  The
  mean-path flux is an imposed relaxed 1-D profile, not a local diffusion, so the limiter
  does not belong there (being removed; R2d).  Deep, the resolved flow at i = 10 carries
  5.6 F_req at 1.4e6 in R2b -- a plume transient still to be understood (R2a comparison).

## Decision 2026-09-09 ~05:00 (user): MLT is OFF the production path

The luminosity deficit without MLT (0.4-0.6 L at 8x8, ~0.33 at 32x32) is a resolution
effect and never lowered the timestep (production dt 30 s flat for 1e7 s; R0f 58.9 s flat);
the photosphere is at the right T in every no-MLT run.  Everything that DID collapse dt was
an addition: the MLT top handover, WB in the ambient medium, and (suspect) ck_nquad = 2.
Production config `prod15/`: seed cut (vpert_rmin 1.6e12), ambient medium + wb_rmax,
extended opacity table, direct RT source, no clipping, ck_nquad 1, mlt_alpha 0.  Gated on
R3_full15 (job 194483): the same at 8x8 to 1.5e6, to see the background's free fall arrive
at the star (~4e5 s) -- never yet run past 3000 s.  MLT stays in the code as a possible
deep-only accelerator (it carried L cleanly from 0.15 to 0.65 R).

## No seed at all? (R4_noseed, job 194486, launched ~05:10)

User's suggestion: vpert = 0 and let truncation error be amplified by the instability.
Expected: surface e-fold ~8e4 s (sup 0.1-0.5, g 16, H_p 1.4e10) -> ~25 e-folds from
~1e-10 imbalance = ~2e6 s to saturation (vs ~1e6 seeded); deep e-fold ~4e6 s -> ~1e8 s,
which is the physical spin-up the seed cut already imposes.  The RT cooling at the top is
a genuine forcing, so the surface should start on its own.  Advantages: no per-block
sin(3x) sin(2xi) sin(eta) imprint, no injected entropy advection, "surface first" arises
from the physics.  Risk: grid fingerprint (seams/vertices) in the first cells.  R4 = R3
config with vpert 0; agent compares onset time, seam fingerprint, deep quiet.  If the
surface reaches v_rms/v_mlt >= 0.3 by ~2e6 s, production goes with vpert = 0.
