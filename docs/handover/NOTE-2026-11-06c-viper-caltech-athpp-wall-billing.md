# NOTE 2026-11-06 viper -> Caltech: Athena++ wall billing per step (item 3 of TASK-2026-11-06b; CORRECTS NOTE-2026-11-06b)

**Answer: billing is not what makes the cold wall cell, and it does not sit mainly at the wall.** In the wall cell (cell 0) it removes about 7-40 L_* (0.7-4.4e42 erg per ~200 s step, 1-6e-4 of the cell's gas + radiation energy per step), against 21-127 L_* of radiation arriving there. The exception is one A_af1 interval (1-2 d) where it is about zero.

**Correction to NOTE-2026-11-06b.** That note said the billing loss is concentrated in the first three cells and inferred that it could cause the film. The first claim came from L_Bd1, a closed-wall run, and does not hold for your two arms. Please drop both. The cold cell also has little billing, and the audit already shows the same wall pattern with billing removed (below).

## Method
Athena++ e47fabb7 + audit ledger patch, rebuilt on viper; the only extra edit is the times at which per-cell profiles are written. Files are in `docs/handover/athenapp-bsg/wall_1106/`.
- Ledger on does not change the run: all dumps are byte-identical to the double-precision reruns in NOTE-2026-11-06b (A_af1 21/21, A_tc4 16/16).
- Per-cell closure is exact at the printed digits (dE = terms, cells 0, 1, 2 and whole domain).
- Reproduces the audit's L_tc4 at 5 d: billing -11.97e45 (audit -12.0), cells 0-10 -11.36e45 (-11.4), boundary -25.85e45 (-25.8).
- B = billing loss = radT_bill - radT_cons (face fluxes that do not telescope). C = radT_cons, transport with all faces from the final intensities. H = hydro flux + point-mass source. Full sphere; L_* = 5.58e38 erg/s.
- Profiles exist to 7 d (A_af1) and 5 d (A_tc4, NaN at 5.07 d). The 10 d profile was not written. Whole-domain and cells 0-10 totals for A_af1 go to 9.9 d from `ledger.txt`.

## Wall cell (cell 0, r = 20.08), mean per interval
| t [d] | A_af1 B [L_*] | C | B per step [erg] | fraction of cell energy per step | A_tc4 B [L_*] | C | B per step [erg] | fraction |
|---|---|---|---|---|---|---|---|---|
| 0-0.5 | -24.7 | +56 | -2.9e42 | -4.4e-4 | -39.7 | +127 | -4.4e42 | -5.7e-4 |
| 0.5-1 | -18.3 | +29 | -2.1e42 | -3.2e-4 | -17.8 | +67 | -1.9e42 | -2.3e-4 |
| 1-2 | +0.3 | +9 | +3.7e40 | +5e-6 | -7.2 | +21 | -7.7e41 | -8.9e-5 |
| 2-3 | -8.4 | +33 | -9.3e41 | -1.2e-4 | -6.9 | +21 | -7.4e41 | -8.2e-5 |
| 3-5 | -12.9 | +29 | -1.4e42 | -1.6e-4 | -11.5 | +25 | -1.2e42 | -1.2e-4 |
| 5-7 | -9.5 | +34 | -1.0e42 | -1.1e-4 | dead | | | |
Cumulative in cell 0: A_af1 -2.67e45 at 5 d, -3.58e45 at 7 d; A_tc4 -3.19e45 at 5 d, against a cell energy of 6.5-10e45 erg. The cell's gas + radiation energy still rises in every interval, so it stays hot despite the billing.

## The cold cell (cell 1) is not billed cold
In L_*, mean per interval:
| t [d] | A_af1 B | C | H | A_tc4 B | C | H |
|---|---|---|---|---|---|---|
| 0-0.5 | +1.7 | -8 | +7 | -0.4 | -60 | +54 |
| 0.5-1 | +1.9 | -6 | +3 | -12.5 | -26 | +61 |
| 1-2 | +4.3 | -33 | +25 | -1.2 | -67 | +77 |
| 2-3 | +7.2 | -83 | +86 | -12.5 | -44 | +64 |
| 3-5 | +4.5 | -36 | +31 | -3.1 | -7 | +11 |
In A_af1 billing adds energy to this cell until 5 d. The large terms are transport out and hydro flux in, which nearly cancel. At 3 d rho x T relative to the IC is 0.96 (A_af1) and 0.86 (A_tc4): the cell is dense and cold at roughly the original gas pressure.

## Where billing does sit, cumulative to 5 d [1e45 erg]
| arm | total dE | B total | B cells 0-2 | B cells 3-10 | B cells > 10 | C (boundary) |
|---|---|---|---|---|---|---|
| A_af1 | -32.49 | -13.10 | -3.86 (29 %) | -7.90 (60 %) | -1.35 (10 %) | -19.50 |
| A_tc4 | -37.69 | -11.97 | -5.50 (46 %) | -5.86 (49 %) | -0.62 (5 %) | -25.85 |
A_af1 at 9.9 d: B -26.6, of it -24.0 in cells 0-10; C -44.5. So billing is 30-40 % of the energy loss and the boundary term is the rest.
Whole-domain B rate [L_*] per interval from 0 to 5 d: A_af1 -19, +5, +11, -50, -113; A_tc4 -63, -30, -10, -32, -80. Most of the late rise is in cells 3-10, which is also where the disturbance has spread.

## Billing removed: what the audit already has (`athenapp-bsg/fix_1016/NOTE.md`, checkpoint 7; not rerun here)
- C_Bd1 and F_Bd1 (converged, zero billing, closed wall, af1) develop **the same pattern**: cell 0 hot, a cold dense cell above it (the audit names cell 2, T 1.3e6 -> 5e5 K, rho x2 at 2 d). F_Bd1 dies at 2.7 d. The audit's verdict: "not caused by the billing or the criteria".
- F_tc4 (his BC, af0, converged) lives past 22 d, where the default dies at 5.07 d. dE -31.1e45 at 5 d, all at the boundary. So billing does matter for how long the run survives in that arm.
- F_af1 (his BC, af1, converged) goes the other way: dt collapses at 8.97 d, with +10.9e45 of energy resets, while default A_af1 runs on.
These are the closest evidence, but not your A_af1/A_tc4 profile question: their wall N^2 was not computed here, and the C/F_Bd1 wall is closed, not his BC.

## What I would conclude
- Your first branch ("the difference to chase is where the cooling lands, i.e. billing") is not supported. The wall pattern forms with billing absent in the closed-wall runs, and cell 1 in your two arms is not cooled by it.
- The wall instability looks like a property of his discretisation at the thick wall, in line with NOTE-2026-11-05e. It does not appear to need the billing.
- This is all 1-D, so why his 3-D stays quiet remains open. A seeded multi-D Athena++ run is still the direct test, and it can be run on request.
