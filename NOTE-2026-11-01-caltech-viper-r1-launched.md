# NOTE 2026-11-01 Caltech -> viper (USER): round-1 WASP-121b LAUNCHED on Caltech; please push the other four planet packages as soon as ready

## Running on Caltech (alloc 5719508, hpc-sm-01-05, 1 H200 each)
- **r1m (10 G resistive MHD + tide, 1x):** started 00:18, deck w121r1_1101/mhd/r1m.athinput.
  - Your s10fresh keys on your hfresh_1x setup plus the Roche anchor (ap = x1min 1.144050e10, grav 1132.4682) and your Roche grid (x1max 2.232001e10).
  - <mhd>/am_vector = true + problem/am_transport = false.
  - Uses the explicit tangential tide until roche_wb_tan supports MHD (being built).
  - Cost on STS: ~178 ms/cycle and ~89 min/rot, so ~18+ days for 300 rot.
  - The resistive limit sits at 1e-3..3e-2 bar on the nightside (eta at the 1e15 cap in the thinnest radial cells), so we are scoping implicit Ohmic resistivity.
- **r1h (hydro control + tide, 1x):** started 00:19, with problem/roche_wb_tan = true. The well-balanced tangential tide is on fin51-r-fix 09487e35; hydro only for now. ~5 min/rot.
- Both still have the old hydrostatic top. The one-way outflow top (topbc-1101) switches in at a restart once gated.

## User: hydro runs of the OTHER planets also go on Caltech, on free GPUs (1 per run, 4 per node)
Please push WASP-76b, WASP-178b, TOI-1518b and WASP-189b as soon as each is ready. Per planet:
- 1x deck on the hfresh defaults;
- inputs with md5;
- Roche anchor inputs (M_p, R_p, p_ref, a, P, Teff, R_*, M_*) or your roche_anchor.py output;
- polar lobe radius and expected dayside r(1e-6 bar) with the tide (see NOTE 61124aaa);
- the UV pseudo-band verdict for the A-star hosts.

We build the tide decks (roche + roche_wb_tan) here with the same tool.
