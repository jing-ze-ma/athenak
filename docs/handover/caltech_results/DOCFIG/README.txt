Palette sheets (NOTE 2026-11-26o item 1), one snapshot: bsg final dump 00114, t = 57.0 d; 2 x 2 panels, 1200 px wide, plot_style.
palette_temperature.png, palette_density.png: cmasher ember, cmasher sunburst, cmasher cosmic, cmocean thermal (theta-r slice, phi = 30 deg, 20-62 Rsun)
palette_vr.png: cmocean balance, cmasher prinsenvlag, cmasher fusion, crameri vik (shell r = 40.70 Rsun, limits symmetric +-64.8 km/s)
palette_flux.png: crameri lajolla, crameri batlowW, and two added to fill the sheet: cmocean solar, cmasher amber (comoving F_r at each column's tau_R = 1 over its mean; 1-99 percentile limits 0.41-1.61)
Versions: cmasher 1.9.4, cmocean 3.0.3, cmcrameri 1.8 (own venv). Scripts data.py (extract), sheets.py (plot).

Round 2 (USER's reply, 2026-11-25 ~12:15): r2_slice_density.png, r2_slice_temperature.png (left (a): delta/<.> from the shell mean in cmasher fusion, cmocean balance, crameri vik; right (b): log zoomed on the convection zone in cmocean thermal, cmasher ember), r2_lateral_density.png, r2_lateral_temperature.png (shells 40.70 and 43.52 Rsun), r2_flux.png (Blues, Blues_r, cmocean dense, cmocean ice; same 2-98 percentile limits 0.473-1.500), r2_vr_fusion.png (settled). Floor gas (density <= 1.5 x floor) grey. Limits (percentiles of the convection zone 38.5-44.0 Rsun): r2_limits.json. Scripts r2data.py, sheets2.py.
