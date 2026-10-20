# NOTE 2026-10-20 caltech -> raven: 3-D accretor page from the fresh tracer run acc3d_trc_1019 (TASK-2026-10-20-caltech-acc3d-page)

URL: https://claude.ai/artifact/U37ZLacbA1JBKRaFz2MzQa (updated in place, title "RY Per Accretor 3-D").

## Build
- Bundle acc3d-page-1020 (bsg-files-1009 @ 4128e6b4), run unchanged with env only, on the Caltech LOGIN node (the expansion CPU queue
  is ~4 weeks), NP 4, INLINE=1. RUN = /resnick/groups/carnegie_poc/jingze/acc3d_trc_1019 (read only, live link 5088058),
  OUT = /resnick/groups/carnegie_poc/jingze/acc3d_trc_1019_page.
- EOS table: existing Caltech copy accretor_1008/rhd/eosdump/eos_x0.70_y0.28_logd-18.txt, md5 e548d9c2777b8b63f979c024ab2dd564
  (= Raven's; no GPU regeneration needed). Rosseland table md5 d81b7a46 (bundle default = run files/ copy).
- ffmpeg: spack ffmpeg-6.0 (libx264). Python: python/3.11.6 + py-matplotlib/3.8.2 modules.
- 3-D: ParaView 5.11.2 pvbatch (osmesa module) renders NOISE on the login node ("GLEW could not be initialized"; one worker also
  SIGABRTs in vtkOpenGLTexture). Replaced by a venv with `pip install vtk-osmesa` 9.3.1 (wheels.vtk.org):
  /resnick/groups/carnegie_poc/jingze/venv_vtk/pvpy.sh (wrapper unsets PYTHONPATH, else the module numpy leaks in). Frames correct.
- Runtime: first full pass 203 s (mk_prod 15 s, mk_tracer 100 s, 2-D frames 45 s, reduce3d 16 s); 3-D re-render with the venv
  + mp4 + page 28 s. Refreshes: OUT/refresh.sh (incremental, same env, login-safe NP 4, prints a 6-line summary + TEXTCHECK).

## Coverage
- hst to t 0.1092 (cycle 20100, 0.231 P_orb; stream on since t 0.0502, 0.48 d); dumps: 22 hydro_w (t 0 - 0.105), 11 m1.
- Frames: 2-D movie 11 (paired hydro+m1 dumps, every 0.01), 3-D movie 22.
- publish/: index.html 1.82 MB; movie/acc3d_eq.mp4 90.7 kB, acc3d_eq_poster.png 191 kB, acc3d_3d.mp4 123 kB, acc3d_3d_poster.png 209 kB.

## By-eye checks (sect. 4)
- Header: t, phase, Caltech layout, binary prodmerge-1019b @ 13206541 (H200 nofma 201918b0), deck v11amwtr (9d1b56bb), chain,
  queue (5088058 running, 59-63 pending) all describe this run.
- Run checks: NaN 0, FATAL 0, NON-CONVERGED 0, hst NaN 0; 11 dt retries (t 0.05 stream switch-on, t 0.071-0.072 contact);
  dt 1.4e-5 static -> 2.9e-6 after contact (t >~ 0.072), flat since. Picard mean ~3.2 (max 8).
- Mass budget: residual |r| <= 8.3e-10 (round-off; Mdom 1.14e6, Min 15.10). Mdot_in 254.92 (nominal 254.84); MR at r_meas now 164.
- Energy: the page's residual (dE - Ebnd - Erel + Emsk) = 2.43e5 = 3.5e-2 of |Ebnd| 6.85e6 -- ABOVE the 1e-3 gate as plotted.
  Without the +Emsk term it is 1.46e4 = 2.1e-3 of |Ebnd|, dominated by the one-time +2.1e4 offset at the first output (see
  NOTE-2026-10-20-caltech-acc3d-trc-budgets: Emsk is already inside Ebnd in this deck). So the bundle's mk_prod.py line 90 probably
  double-counts Emsk; please confirm which convention v11amwtr's hst uses. Not fixed here.
- Lz: dJdom - Jin = 5.19e4 (dJ 1.03e5, Jin 5.12e4); the static-phase linear drift (~7e5 per code time) continues into the stream
  phase, and the residual turns over after t ~0.09 (5.8e4 -> 5.2e4). Same open question as the budgets NOTE.
- theta walls: wall flux fraction 9.8e-7 (closed).
- Tracers: dNz0/dt = dMs0/dt = 254.92 (nominal 254.84); s0 budget residual 1.4e-13 (round-off), s1 1e-21; s1 mass 1.3e-7 (negligible).
- s0-weighted entropy of the accreted gas (r < 10.3 R_sun), s_gas [k_B/baryon], mean (p10 / p50 / p90):
  t 0.075 (+0.20 d): 39.8 (36.3 / 40.5 / 41.9), M 0.71
  t 0.090 (+0.32 d): 38.3 (36.9 / 37.5 / 41.3), M 4.35
  t 0.105 (+0.44 d): 38.2 (37.1 / 37.9 / 40.6), M 8.18 of 14.44 s0 in the domain
  Stream gas in flight: s 39.1. Star at t 0.05 (|lat| < 10 deg): s 38.6 at tau 100 (r 9.17), 40.3 at tau 2/3 (r 9.36).
  => the tracer value 38.2 (p10-90 37.1-40.6) agrees with run11's top-dM estimate 38.5 (38.0-39.1), with a wider high tail.
- Where the s0 gas sits: a thin equatorial skin at the accretor surface, r p10/p50/p90 = 8.69 / 8.88 / 9.33 R_sun (R_acc 9.0),
  i.e. it has sunk to ~0.1 R_sun BELOW the original surface; log10 tau_R p10/p50/p90 = 2.7 / 3.6 / 3.9 (buried at tau ~4e3), confined
  to |z| <~ 1 R_sun near the impact meridian and spreading in phi; dM/dr peaks at r ~8.8. The rest of the domain s0 (~6.3) is the
  stream / Roche-lobe fill at r 9.5-14.6 (s0 = 1 there).
- Omega(r): accreted layer spun up to Omega/Omega_orb ~1.2-1.6 at r 8.7-9.2 (s0-weighted latest up to 2.7 at r 8.8), outer
  sub-Keplerian disc-like flow 2.6 -> 1.2 for r 9.5 -> 14.6; star below 8.5 still corotating.
- Maps / 2-D movie: log rho, log T_gas, log P_rad/P_gas, v_phi as specified; ambient saturates at P_rad/P_gas 1e5 (red), envelope
  ~0.1-10 (blue/white), stream + splash visible at the impact. Cosmetic known issue: bottom-row "(phi N)" title suffix clipped.
- 3-D frames checked: t 0.010 (static, symmetric star), 0.070 (stream reaching the star), 0.105 (hot impact spot, T up to ~10^5).

PUBLISHED 10-20 07:15 PDT to https://claude.ai/artifact/U37ZLacbA1JBKRaFz2MzQa (version 13; index.html 1.82 MB + 4 movie files replaced run11's).
