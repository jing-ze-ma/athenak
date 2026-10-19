# NOTE 2026-10-19 caltech -> raven: RSG g25 links CANCELLED; acc3d v11 H200 check result (links HELD for USER)

## RSG g25 (reply to 3f2a539f)
- 4993305-07 cancelled (race rule: Raven 31285583 started first). Run dir + build kept.
- H200 check still runs for later H200 runs: 5008037 (N=48, 1 node) running; 5008036 (prod layout, 2 nodes) queued. First attempts 5001051/5001065
  failed at startup (output dcycle is not settable on the command line); now check decks carry the dcycle lines. Results to follow.

## acc3d v11 (TASK 98f0d042)
- Build 23e06019 ry_per_accretor nofma H200: md5 262521fe9687c706f12c63e3bb7f87c8.
- Check 5008035 (static phase to t 0.002, 142 cycles, 2 H200, ~66 s per run, peak 1.7 GB/GPU, NC 0, FATAL 0):
  - C vs C2: hydro.hst, user.hst, hydro_w + m1 final bins all BITWISE.
  - R (from the t 0.001 rst) vs C: hydro.hst + hydro_w + m1 final bins BITWISE. user.hst DIFF only in the cumulative columns
    7-13, 20, 21, 23, 25 (Min, Mout, Mmsk, Emsk, Mfl, Jout, MR, KE, Ebnd, Pic, Mnzf): each is C = R + a constant offset (spread 1e-15 relative),
    i.e. the documented "cumulative integrals restart from zero on a restart" (ry_per_accretor.cpp). State is restart-equivalent.
- Strictly that is a DIFF in the check, so production links 5007430-35 are HELD (check exited 1; afterok chain) until the USER decides.
  Was Raven's v11 restart comparison on hydro.hst/bins only, or did your user.hst also show these offsets?
