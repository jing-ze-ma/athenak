# Plot style (USER 11-24): the AREPO-IDORT Zenodo notebooks' style, for ALL plots on every machine

Files: `mesa2.mplstyle` (md5 859b82680fd5b16044da6d51803bd317, as shipped with every notebook of
Zenodo 15032138) and `idort_style.py` (`idort_style.use()`).
How the notebooks use it (shock_1d, shadow_2d, wave_2d, diffusion_2d):
- `mpl.style.use('default')`, then `plt.style.use('mesa2.mplstyle')`; text, axes, ticks black; white
  background; ticks inward on all four sides with minor ticks; no legend frame; 300 dpi.
- Figure sizes: 9 x 9 (shock panels), 10 x 10 (slices), 10 x 13 or 10 x 14 (stacked panels, `hspace` 0-0.01),
  12 x 7 (two panels, `wspace` 0.4).
- Labels: axis labels 25-30 pt, tick labels 25 pt, legends 20 pt (wave_2d: `font.size` 20; shock_1d: the style's own).
- Lines: `'k'`, `'r'`, `'b'` (then `'w'`, `'aqua'`, `'lime'` on dark), else the default cycle `C0, C1, ...`;
  reference / exact curves dashed (`'--'`, `'r'` or grey `'0.5'`) or dotted (`':'`, black); simulation points `'k'` or
  coloured `'x'` / `'*'` markers without fill.
- Guide lines for order: dashed `'0.5'` with a rotated `$\propto N^{-2}$` label of the same colour.
- Maps: `'jet'`, colour bar labelled in `\rm`, e.g. `$\rm \log_{10} \rho$`, with tick labels of 25 pt.
- Math labels in `\rm`: `$\rm \tau_a$`, `$\rm E_r$`.
