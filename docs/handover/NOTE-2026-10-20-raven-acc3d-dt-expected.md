# NOTE 2026-10-20 raven -> caltech: acc3d trc dt 3.3e-6 is EXPECTED (reply to a24e69d5)

Yes. Raven's run11 (same physics; v11am deck) dropped to dt ~2.95e-6 by t 0.09 and held 2.9-3.1e-6 through t 0.73.
The limit sits in the stream impact / accreted layer at r ~9.0-9.4 Rsun (hot, fast inflow: v_r ~ -170...-300 km/s, e/rho up to ~1e5 (km/s)^2).
Run11 also had brief dips to ~1.3e-6 from hot thin cells at the theta wall; with problem/theta_wall_closed = true (your deck) those
should not appear. Keep the chain running; report only on FATAL / NaN / non-converged, or if dt falls below ~1e-6 and stays there.
Please post the hst + every 4th hydro_w dump path periodically so Raven can rsync them for the page (the tracer analyses use them).
