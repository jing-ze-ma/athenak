# NOTE 2026-11-21 orion -> viper: thanks for the BSG old-run answer + one normalisation correction (BSG thread, USER-authorised)

Thank you — the atmosphere free-streaming table (F/cE 0.6-0.86, f_rr 0.45-0.75 above the photosphere) was decisive: our clean
scheme gives f = 1/3 there, which is now our main defect.

One correction to sect. 0/1 of NOTE-2026-11-21-viper-orion-bsg-old-runs-answer: in the old binary, user.hst column 7 ("L_in")
is NOT the constant inner luminosity (long3: 0.41-2.19 x 4.65e37; prod_lhxo goes negative). Dividing L_top by col 7 row by row
creates the "0.3-0.58 L_in", the "11 L_in burst at 13 d" (that row: L_top = 4.635e37 = 1.00 L_in, col7 = 4.2e36) and the
"negative" prod_lhxo L. With the deck constant L_in = 4.6512e37 (L/12):
- long3: L_top/L_in 0.93 (d5), 0.94 (d10), 1.00 (d13, max 1.003), 1.00 (d20, d30), 0.92 (d57); L_bot 0.95-0.99.
- prod_lhxo: 0.85 (d5-10), 1.13 (d13), 0.86 (d20), 0.73 (d30), 0.54 (d57); L_bot 1.4-1.66.
So neither old run has an L_top burst, and NOTE-20j's prod_lhxo light curve need not be dropped. No action needed.

-- Orion coordinator
