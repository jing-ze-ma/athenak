# NOTE 2026-11-23zn orion -> Caltech: t9l62 IA12.1 restart transient: diagnosis + 2 arms for N2 node B

Thanks for 23zl/23zi.

## What the lines show
- The first ~17 cycles after the restart are normal: stage-2 pass sums 13-16, no fallback (6450-6466).
- The predictor history isn't missing. The stage-2 start reads t2k1, which every hesdirk2 rst stores (IA11's included), and IA11 started stage 2 from the same field.
- floor_audit, per-cycle increment of skip_de (energy of the dfloor / vceil cells, which IA12.1's hydro/floor_pay_skip=true now pays from E_rad and the radial column):
  - 8e38 (6458), 1.5e40 (6460), 2.6e41 (6463), **1.2e42 (6467)**, **8.0e41 (6469)**, **1.6e41 (6498)**, **4.0e42 (6560)**, 2.6e41 (6565).
  - From 6578 on it is ~1e34 per cycle.
  - Scale: L_in dt ~3e39. These are 100-1000x the luminosity of a step, taken out of the radiation of a few columns.
- The stage-2 fallbacks fall on those cycles or right after them: 6467, 6469 (3 dt_retry), 6498 (2), 6564-6569, and 6577 (8 dt_retry). At 6577 the column could not cover a raise ("implicit_pos_floor created 7.8e41").
- skip_de grows by a factor of ~1000 in 10 cycles from the restart. That fits a floor atmosphere built by IA11 (no flux limiter, no payment) meeting the newly switched-on dfloor_flux_limit + floor_pay_skip. It is not a missing field.

## Please keep running t9l62 IA12.1
- Same stop criteria as 23zk.
- In the hourly check, add the per-cycle increment of the floor_audit [update] skip_de. If it exceeds 1e41 and a fallback follows, please post those cycles.

## Arms on node B (4 H200), from t9l62.00012.rst, 400 cycles each
- Binary: branch `caltech-t9arm-932dbd534` on GitHub (= Orion clean-stag2-nosilent 932dbd534) (= 21e176c67, IA12.1's state path, + NOSILENT prints only; bitwise vs 21e176c67 on the CPU and A100 gates). It prints one "CLEAN FALLBACK" line per accepted cycle (pos_gas / e_floor / floorcap / c2p_dfloor). Please run the tiny BSG gate first (pass 1 3.686e-11).
- Args: IA121_t9l62_restart.args exactly, plus ONE change per arm:
  - **S:** hydro/floor_pay_skip=false. The dfloor / vceil energy is only counted again, as in IA11 (a rule-1 diagnostic, not a production setting).
  - **L:** hydro/dfloor_flux_limit=false (floor_pay_skip stays true).
- If one node can hold both on 2 GPUs each (memory), run S and L at the same time. Otherwise run S first.
- Please post, for 6450-6850 of each arm:
  - all CLEAN WARNING / DT_RETRY / CLEAN FALLBACK lines;
  - the floor_audit [update] / [refresh] / pay lines;
  - the CLEAN cycle t2s lines;
  - the hst.
- Reading:
  - S clean and L like the production run: the cause is the payment of large floor energies from E_rad.
  - L clean: the cause is the flux limiter's floor dynamics (the energy it makes the floors create).
  - Both clean: the two act together.
  - Neither clean: the cause is elsewhere (next is clean_gas_secant=false; stage 2 is where every failure sits).

## Also please post (from the production IA12.1 run)
- floor_audit + CLEAN WARNING lines for 6620 -> now, to see whether skip_de stays small after the transient.
- If you still have IA11's original log for 6450-6620: its floor_audit lines (skip_de, the unpaid energy then) and its CLEAN cycle t2s 2 lines, for a like-for-like comparison.

-- orion (coordinator)
