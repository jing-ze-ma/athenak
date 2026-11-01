# Roche top check for every dhj_planets_1101 package (11-01, after Caltech NOTE 61124aaaf)

Script `roche_top_table.py` -> `roche_top_table.out` (no runs). It uses the tide-on roche_anchor values of each deck
(ap = x1min, grav, omega, roche_a), and the potential and IC column classes of roche_anchor.py c8aac483c.
- **r_polar** is the polar lobe radius, where Phi_tot on the z axis equals Phi(L1).
- **r(1e-6) t=0** is the 1e-6 bar level on the tide-on equipotential of the IC column.
  - Check: the WASP-76b CPU fresh start with the tide on (gate/run_x224) gives 1.7367e10 substellar and 1.5201e10
    at the pole; the tool gives 1.7379e10 / 1.5199e10.
  - For W121 1x the tool gives 1.637e10 substellar; Caltech quotes 1.564e10. The 3x anchor and IC give 1.597e10.
    That difference is unexplained.
- **Evolved estimates** carry Caltech's W121 rise above their t=0 isobar (+0.166e10 after 1-3 rot, +0.236e10 evolved),
  scaled by each planet's substellar p_ref..1e-6 bar height. Starting from our higher t=0 value makes this conservative.

| planet | x1max [cm] | r_L1 | r_polar | x1max/r_L1 | **x1max/r_polar** | r(1e-6) pole t=0 | **r(1e-6) substellar t=0** | 1-3 rot est. | **evolved est.** | x1max/evolved | margin |
|---|---|---|---|---|---|---|---|---|---|---|---|
| WASP-121b 1x | deck 2.1504e10 | 2.4225e10 | 1.5595e10 | 0.888 | **1.379** | 1.384e10 | **1.637e10** | 1.803e10 | **1.873e10** | 1.148 | +2.8e9 |
| WASP-121b 1x | roche grid 2.2320e10 | 2.4225e10 | 1.5595e10 | 0.921 | **1.431** | 1.384e10 | 1.637e10 | 1.803e10 | 1.873e10 | 1.192 | +3.6e9 |
| WASP-76b 1x | **FINAL deck 2.4000e10 (tide-on refit)** | 2.8043e10 | 1.8034e10 | 0.856 | **1.331** | 1.520e10 | **1.738e10** | 1.933e10 | **2.016e10** | 1.191 | +3.8e9 |
| WASP-76b 1x | 2.2400e10 (superseded) | 2.8043e10 | 1.8034e10 | 0.799 | 1.242 | 1.520e10 | 1.738e10 | 1.933e10 | 2.016e10 | 1.111 | +2.2e9 |
| WASP-178b 3x | deck 1.8985e10 | 5.1795e10 | 3.3340e10 | 0.367 | **0.569** | 1.395e10 | **1.413e10** | 1.497e10 | **1.532e10** | 1.239 | +3.7e9 |
| WASP-178b 1x | deck 1.9325e10 | 5.1795e10 | 3.3340e10 | 0.373 | **0.580** | 1.408e10 | 1.427e10 | 1.519e10 | 1.558e10 | 1.240 | +3.8e9 |
| TOI-1518b 1x | deck 2.0078e10 | 3.6682e10 | 2.3626e10 | 0.547 | **0.850** | 1.491e10 | **1.562e10** | 1.673e10 | **1.720e10** | 1.168 | +2.9e9 |
| TOI-1518b 1x | roche grid 2.0309e10 | 3.6682e10 | 2.3626e10 | 0.554 | **0.860** | 1.491e10 | 1.562e10 | 1.673e10 | 1.720e10 | 1.181 | +3.1e9 |
| WASP-189b 1x | deck 1.5953e10 | 4.7963e10 | 3.0879e10 | 0.333 | **0.517** | 1.270e10 | **1.286e10** | 1.359e10 | **1.390e10** | 1.148 | +2.1e9 |
| WASP-189b 3x | deck 1.5673e10 | 4.7963e10 | 3.0879e10 | 0.327 | **0.508** | 1.256e10 | 1.271e10 | 1.334e10 | 1.361e10 | 1.151 | +2.1e9 |

r_polar / r_L1 = 0.64 for every planet. For W121 our 1.5595e10 matches Caltech's 1.5596e10.

## Recommendations (x1max)
- **WASP-121b 1x: no change.** The deck (2.15e10) and the roche grid (2.23e10) are on open polar equipotentials; Caltech's
  outflow1 top is meant for exactly this. The dayside margin (x1max/evolved 1.15-1.19) is the reference Caltech accepts.
- **WASP-76b: DONE (user decision 11-01).** x1max is back at 2.40e10, refitted on the tide-on anchor (`wasp76b/grid/roche_x240`,
  max dr/target 0.984), with roche = true. It is gated (CPU c8aac483: r(p_ref) +1.6e-4 / -2.7e-4; GPU apudev 12529652).
  The 2.24e10 deck is superseded. The code's t=0 substellar 1e-6 bar level is 1.7364e10.
- **WASP-178b 3x and 1x, TOI-1518b 1x (deck and roche grid), WASP-189b 1x and 3x: no change.**
  - The whole top shell is inside the polar lobe (x1max/r_polar 0.51-0.86), so outflow1 is not needed.
  - The dayside 1e-6 bar margins are 1.15-1.24. Caltech's W121 1.49e10 failure was a top 0.95 r_polar below the evolved
    isobar. These decks keep the isobar inside, with the same or more margin than W121.
  - The W189 margin (1.15, +2.1e9) equals W121's in ratio, and its scale height is smaller (g 2100).
