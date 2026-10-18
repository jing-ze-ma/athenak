# NOTE 2026-10-18 raven -> caltech: static-phase mid-latitude infall -- SAME on Raven v10; continue

Raven v10 static (runs/ctc10/static, 4 A100, same binary/deck) dt_report: dt flat 1.40e-5 to t 0.032, then the thin-atmosphere
cells take over: t 0.0338 dt 9.9e-6 (r 9.41, v_r +129 km/s), 0.0358 8.5e-6 (r 9.59, +158), 0.0378 1.12e-5 (r 9.11, -109),
0.0398 9.8e-6 (r 9.13, -127), 0.0418 8.2e-6 (r 9.14, -157), theta ~108 deg = off-equator, both signs = the thin-atmosphere
breathing/von Zeipel circulation (rho < 1e-6 code), not a deep-envelope problem. The OLD barotropic IC (ctc9) shows the same
phenomenon slightly later and milder (from t 0.038, -94..-128 km/s near the equator, dt 1.0-1.2e-5). So v10 costs ~1.5-2x in dt
in the static phase; no retries/NC anywhere. Continue (USER-chosen IC); the stream switches on at t 0.05 as planned.
