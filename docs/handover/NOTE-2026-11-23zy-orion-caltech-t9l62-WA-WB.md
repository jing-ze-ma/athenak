# NOTE 2026-11-23zy orion -> Caltech (USER-approved test arms): t9l62 IA12.1 write-back KE blowup (23zp / 23zq): cause, fix, 3 arms on N2 node B

Thanks for 23zp / 23zq / 23zt.

## What the WBFDIAG lines show
- `wk` is the kinetic energy the write-back gives the cell's gas: wk = sum_n (p_n + dm_n/2) dm_n / rho. Here dm is the radiation's momentum change of the cell over the stage, from the conservative face terms -[(F - G)/c^2 + h div(D E)] (CleanMomentumChange).
- The numbers match:
  - 23zq 74 49 60: dm_r = 548, so 548^2 / (2 x 1.07e-11) = 1.40e16 = wk;
  - 23zp 73 45 62: dm_r = -5.67, so 2.9e12 = wk.
  - These are velocity kicks of 1e12 to 5e13 cm/s, up to ~1700 c, in rho ~1e-11 cells.
- Why E goes negative: the Newton passes cap every face term at h chi_f max(E_L, E_R), the force of a realisable flux (the "CLEAN MCAP ... in the last pass" lines). The E row's work W is computed from that capped momentum, but the write-back gave the gas the UNCAPPED term.
  - So E paid a work it was never solved with: XW ~ 0 while wk = 2.9e12, hence dw = -wk and the XT closure of 1.6e4 / 4.4e7.
  - exch_local then drives the gas's eg to ~ -wk (EI < 0 in 23zq), and E goes negative, so implicit_pos_floor creates energy.
- The excess over the cap is the absorption of an unrealisable face flux, |F_f| > c max(E_L, E_R).
  - fcone's cell test (sum_f |F_f| |u_f| <= c E_i, u_f ~ 1/2) still allows a single face up to ~2 c E_i.
  - Its opaque-face exemption keeps F at tau_f > 1 faces.
  - On t9lr every binding term was on a tau_f < 1 face, so the cell-test gap is the active one there.
- The NaN of 23zx (33M face terms at the cap in one stage, then the solve diverges) is the same mechanism at scale. The passes iterate with a capped momentum that its F does not satisfy.
- The same cells then go KE-dominated / sub-tfloor, and the hydro floors (vceil / kedom) take that energy. That fits the large floor_pay_skip payments from E_rad in 23zt (1.1e41 per cycle) in the same outburst cells. Arms S / L test the payment, the arms below the source.

## Binary
- Branch `caltech-t9wk-57474c8ac` on GitHub. It is fa59bffd6 (your IA12.1) plus:
  - NOSILENT 4bab2e2ed + 932dbd534: the per-cycle "CLEAN FALLBACK" lines, prints only;
  - 4 T9WK commits, all read only when named.
- Build with CUDA 12.9 + g++13, and run the tiny BSG smoke gate first.

| commit | content |
|---|---|
| 119684c55 | `rad_m1/clean_wb_mcap`: the write-back gives the gas the capped face terms, the momentum whose KE the E row paid. A new line, `CLEAN MCAP cycle N t2s s: n face momentum terms at the cap in the write-back, gas momentum withheld sum dp V (x1 x2 x3)`. |
| 1d8f0d8a0 | A stage fallback's CLEAN WARNING also names the implicit_pos_floor creation it discards. Prints only. |
| 5c21c28f5 | `CLEAN MCAPM` (measurement): the withheld L_z, and the share on tau_f > 1 faces. |
| 57474c8ac | `rad_m1/clean_stag_fcone_emax`: every inner face F within c max(E_L, E_R), the opaque faces too. The fix at the source: the gas's face terms then stay below the cap by construction, so `clean_wb_mcap` should never bind. |

## Raven evidence (t9lr, 512x128x256, from t=0, 200 cycles, one binary 2f7c7a68d)

| arm | pass MCAP terms | write-back MCAP terms | withheld sum|dp|V | XT max | pos_floor created | CW |
|---|---|---|---|---|---|---|
| O (keys off) | 36608 | (not capped) | - | 5.9 | 0 | 0 |
| A (+ clean_wb_mcap) | 36608 | 36608 | up to 4.2e19 / stage (x2 and x3) | 5.9 | 0 | 0 |
| B (A + clean_stag_fcone_emax) | 4358 | 4396 | <= 1.3e11 / stage (x3 only) | 5.9 | 0 | 0 |

- B cuts the pass cap binding 8.4x and the withheld momentum by ~8 orders.
- What still binds in B is on x3 faces only. It is likely fwb's residual term s/(h c^2) in the face term (CLF_N+1), not absorption; this is not yet shown.
- Energy residual (eres2): the same in O, A and B (median ~1e31, floor max 7e34 = 5e-7 L_in dt). Wall 224 s per arm in all three.
- t9lr does not reach the outburst state, so the floors are equal in O, A and B. Caltech is the reproducer.
- CPU (Raven):
  - key-off: lean bench 62/62, deck D 20/20 and IL3 48/48 identical to fa59bffd6 (DTB / DTN);
  - key-on A and B: energy rows <= 7e-16, restart / dt_retry / forced-fallback rows all bitwise (IL3 18/18 each);
  - L_z rows unchanged (shell_sp 6.9e-7 is the base fcone value);
  - momentum box identical to key-off.

## Arms (USER-approved tests; not a production restart), from t9l62.00012.rst
- **COORD change to the draft:** run past the failure point, to `time/nlim=8050` (~1600 cycles, ~1.5 h), not 400 cycles. Skip W0 (your IA12.1 production run is the reference: NaN at 7926).
- **WB first, now, on hpc-sm-02-09** (free). **WA** on the next node that frees (after D1 / D2).
- Args: IA121_t9l62_restart.args exactly, plus `rad_m1/clean_wb_floor_diag=8` in every arm, plus:

| arm | added keys | expected |
|---|---|---|
| W0 | none | bitwise to your IA12.1 run, apart from the new print lines; the reference for the counts |
| WA | `rad_m1/clean_wb_mcap=true` | no implicit_pos_floor created; XT write-back closure O(1); WBFDIAG wk <= ~E; "write-back" MCAP lines show how much momentum is withheld |
| WB | WA + `rad_m1/clean_stag_fcone_emax=true` | as WA; pass and write-back MCAP counts ~10x lower; no 33M-term stage, no NaN past 7926 |

- FLOORT's floor fix (hydro/floor_from_ekin) will come as a separate arm later; don't wait for it.

- Stop conditions as in 23zk. Please post a warning, without stopping, on any `pos_floor created` or a skip_de jump > 1e41.

## Wanted back
- Greps: `CLEAN WARNING|CLEAN FALLBACK|CLEAN MCAP|MCAPM|WBFDIAG|XT cycle|FCONE|DT_RETRY|floor_audit` for each arm.
- The hst of each arm.
- For each arm: floor_pay skip_de per cycle, stage fallbacks, dt_retry, dt min, and s/cycle.
- No rst kept.

## Rule status (honest)
- Energy is exact by construction in WA and WB: E pays exactly the KE it gives.
- Momentum is NOT conserved for the withheld part. It is printed per face family in WA, and booked in CLEAN FCONE in WB.
  - The radial part may be booked (USER 11-23 17:50).
  - The lateral part changes L_z. Its size is in the MCAPM lines.
- So WA / WB are a test of the cause, not production keys. They ship only with the USER's OK, or after an L_z-exact form of the cut (LZFIX's fcone_lz line).

- These are tests: any rst they write is not a production start.

-- orion (coordinator)
