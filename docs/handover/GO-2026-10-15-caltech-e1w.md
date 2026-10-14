# GO (USER, via Raven coordinator, 10-15): E1W remapped narrow-wedge production

USER chose the remap start of the 2x-lateral pi/6 wedge E1 (option 2). Gate G0/G1/G2 PASS (53ab04cb) -> submit
e1w_1014/prod/submit_links.sh as planned (1 H200, QoS normal 7-d link + 2 afterany, tlim 1e9, UCX_TLS=^xpmem, binary 6d4215c0).
Reports side by side with e1s at matching times (x3.8638 for wedge totals), incl. Picard mean, BE-step / stage-fallback and DT_RETRY
counts (G2 showed Picard mean ~21, max 152, 46 BE steps in 300 cycles -- the same late-ramp stiffness e1s shows), s/cycle at
production cadence, lateral cells per H_p in the outer CZ. Flag at once: infall < -5 km/s, rsurf > 1.15, DT_RETRY bursts, Picard mean > 40.
e1s untouched.
