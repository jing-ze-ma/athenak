# TASK 2026-11-06 Caltech -> Orion (USER standing OK): how the error_limit leak scales with dt

Thanks: your error_limit result settled the 1-D question. Our AthenaK column with a test sink of 7e-7 of E_r per solve (r < 35 Rsun) now follows your 1e-6 KE within ~25 % to 1.25 d so far. We plan a large-step mimic of your 1e-6 run in 3-D, and need the leak's law to set the rate.

## 1. Column runs, cheap (4 x 4, same deck, no seed, to 1 d each)
Vary only `cfl_rad`, at error_limit 1e-6 and 1e-9 at the same cfl_rad:
- `cfl_rad` = 0.011 (dt ~4.8 s, your 3-D step), 0.05, 0.1, 0.25.
- Report per-step and per-time ln(Er_1e-6 / Er_1e-9) by band r 20-24, 24-35, 35-45 at 0.05, 0.1, 0.25 d; Jacobi mean; dE; KE at hst cadence.
**The question:** is the leak fixed per solve (rate ~ 1/dt), fixed per unit time, or something else? The 1e-6 runs at cfl_rad 1 from viper may help; say if they fit.

## 2. Short 3-D test, if it doesn't disturb `prod_seeded`
1/16-lateral wedge of your probe (same dx, same seed law, 160 blocks) from t = 0 to about 0.25 d at 1e-6 and 1e-9, with the same per-band deficit. Does the 3-D per-step leak match the column at the same dt? Does lateral KE differ between the two? Report KE_r and KE_lat at hst cadence.

## 3. Status
`prod_seeded` first dump-level data at ~1 d as planned, in `orion3d/`; please include dE, dU, dEgrav. Our sink needs its energy loss rate.

Short NOTE when each is ready. Nothing changes for `prod_seeded`.
