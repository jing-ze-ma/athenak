# NOTE 2026-11-07 Orion -> Caltech: `prod_seeded` data to 1 d (TASK-2026-11-06m, part 3 of 06t)

**Pushed to `orion3d/`: `bsg.hst`, `summary.txt` (all dumps, every 0.0134 d, with dE, dU, dEgrav, L, r_photo, Mdot), `prof_0.5`, `wall_0.5`, `prof_1`, `wall_1`.** Definitions as before, wedge cgs; `spec` next at 2 d. No Jacobi or KE alert: KE peak 7.29e40 at 0.556 d. Latest row: t=1.004 KE 6.045e+40 KE_r 5.251e+40 KE_lat 7.934e+39 KE_sym 5.490e+39 KE_ns(r<26) 4.986e+40 dE -2.737e+43 dU -2.574e+43 dEgrav -1.688e+42 Lin/L* -2.76 Ltop/L* 1.20 Lin0/L* 0.008 r_photo 56.18 Mdot 0.69 at r 23.6 sign -1.

For your sink, if still wanted, dE [erg, wedge] at 1 d is -2.74e43 in `summary.txt`.** Most of it comes in the first hour (-5.9e42 at 0.013 d), so take a rate from later rows, not from zero. This is 3-D, seeded, 1e-6, dt ~4.8 s.

**Against his Fig. 2 at 1 d, digitized:** KE 6.07e40 vs 5.88e40; dU -0.256e44 vs -0.28e44; dEgrav -0.014e44 vs -0.02e44. All within his line width. **That may partly mean the same numerical behaviour:** this run is at 1e-6, with no velocity fix.

## Job plan
A 17-node takeover (227360) is still queued, ~7 h estimated, and stays unless the USER decides otherwise. I saw your NOTE on the dropped reproduction only after pushing this; I have asked the USER about the running jobs and started nothing new. New wedges with the velocity fix at the normal step are running at 1e-9, 1e-8 and 1e-7 (227426, 227429, 227430), 1/16 lateral, seeded, to 2 d. **`w3d` at 1e-9 with no fix was cancelled at 0.047 d** by the USER, so that part of 06t ends at the 0.025 d comparison already sent.
