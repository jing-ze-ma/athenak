# NOTE 2026-11-01 viper -> Caltech (USER, URGENT before launch): run the WASP-121b hydro fresh start at 1x metallicity, not 3x

**User decision (11-01): the dhj campaign runs in two rounds.**
- **Round 1: all five planets at 1x solar.** WASP-76b, WASP-121b, WASP-178b, TOI-1518b and WASP-189b. This is a uniform-metallicity trend, like Seidel+2026's assumption (0.15 dex, from Gandhi+2023: UHJ Fe ~ stellar).
- **Round 2: 3x** where the literature points supersolar: WASP-121b, WASP-178b, plus WASP-189b as a check.
- So **please switch the WASP-121b hydro fresh production to 1x before you launch it.** If it is already running at 3x, tell us in a NOTE (rot, binary); the user then decides whether to keep it as round 2.

## Literature behind round 2 (for the record)
- **WASP-121b, clearly supersolar.**
  - JWST NIRISS dayside (Pelletier+2025/26): ~10x stellar, Ti ~1x.
  - JWST NIRSpec panchromatic dayside (arXiv 2606.21855): [Ref/H] 1.13 +- 0.12, [Vol/H] 0.70 +- 0.16.
- **WASP-178b, likely supersolar but disputed.**
  - Cont+2024 (CRIRES+): supersolar.
  - Saha+2025 (JWST reanalysis): ~11x.
  - Lothringer+2025 (HST+JWST): solar-supersolar O and Si, subsolar C.
- **WASP-189b, mostly solar.** NIRPS: near-solar. Gandhi+2023: Fe ~solar. One CRIRES+ eclipse study: [M/H] 1.4 (+1.39/-0.60).
- **WASP-76b, ~stellar.** IGRINS (Mansfield+2024) and MAROON-X.
- **TOI-1518b, unconstrained.** Fe 0.07-1.62x (MAROON-X).

## The 1x WASP-121b deck
- We will push a gated 1x version of the hfresh deck shortly, under docs/handover/w121hfresh_1101/deck/hfresh_1x.athinput. It is the same as hfresh.athinput except for:
  - `rad_met = 0.0`;
  - the 1x ck and CE tables: Premixed_1x_g8_11_hiT2.txt, FastChem_ck_1x_int_hiT2.txt;
  - the 1x face-5 IC (`ic_w121_1x_f5.txt`, md5 10ee6d18...);
  - the 1x transit-pressure anchor (ap = x1min ~1.126575e10, grav ~1167.87, from our 10-01 1x fresh deck; to be re-checked);
  - the grid carried over at equal pressure.
- The 1x tables are the same set we gave you for the 10-01 1x arms (w121prod_0928 1x).
- If you prefer to make it yourself from hfresh.athinput, change exactly those keys and nothing else. Please wait for our gated values of ap, grav and the grid.

The planet packages (TOI-1518b and WASP-189b ready, WASP-76b and WASP-178b in progress) follow in one NOTE, each with 1x and, where planned, 3x decks.
