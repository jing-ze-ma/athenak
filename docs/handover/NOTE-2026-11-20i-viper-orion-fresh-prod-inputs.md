# NOTE 2026-11-20i viper -> orion (USER): fresh v20.1 inputs + floor tests; USER: run on Caltech, T9L6.2 first

Reply to NOTE-2026-11-20h. Everything is in `docs/handover/fresh_prod_1120/`: read `README.md` there first. Data files
are gzipped and decks use a `FRESH_DIR` path placeholder; checksums are in `MD5SUMS`.

**USER decisions (11-20 PST):**
1. **Ideal gas for every star, following the papers** (an exception to the general-EOS rule). μ, GM and L are each
   checked against the IC; see the README tables.
2. **Orion runs these on Caltech.** Nothing is launched on viper.
3. **T9L6.2 first.** This is the outbursting LBV model of Jiang+18 (IC init_outburst.txt). Run its corona test before
   any long run.

**Contents:**
- `t9l62/`: the our-wall deck plus its IC. `t9l62/corona_test/` holds the floor test: a static hydrostatic corona
  above 95.5 Rsun, dfloor 1e-16, and a new floor sponge with dmax 1e-15. Suggested tlim is ≥ 8e5 s, past the viper
  paper-BC crash at ~189 h (DTDROP_1115.md included).
- `t19l64/`: the deck plus its IC. `t19l64/corona_test/` applies the same recipe: corona above 31.6 Rsun, dfloor 1e-16,
  sponge dmax 1e-15, tlim 1.8e5 s.
  - It is motivated by the viper 11-19 run: from 30 h, near-empty cells (ρ 5e-19-1e-14, T 2-7e8 K) sat at the 3e8 cm/s
    vkick cap at 37-45 Rsun and in evacuated pockets at 23-25 Rsun.
  - Over the same time, a dense outflow at 1000-3000 km/s left the Fe bump (IC Γ up to 1.13).
- `bsg/`: the Ma+2026 reproduction deck with **μ 0.614**, the value the Linhao profile's own HSE gives; the paper and
  the old deck have 0.62.
- `m35mms/`: Will's paper deck plus the ideal-gas IC (μ 0.6215).
- `m35zams/`: the new ideal-gas IC (μ 0.6220, HSE 4e-8, Γ max 0.88). **No deck yet:** take the grid from
  GRID_REFERENCE (S22 T42L5.0, 384x128x128) and the EOS/IC/BC keys from m35mms.athinput.

**OPEN (USER):**
- BSG GM. The deck has 19.78 Msun from his deck, against the profile's M(r) of 19.82-19.99 Msun over the domain;
  19.90 would keep the error within ±0.5 %.

**Please report to viper, for each corona test against its baseline:**
- dt history;
- implicit_wb_vkick_reject counts (kick solves, cells, stage rejects) and dt COLLAPSE / NOT ADMISSIBLE counts;
- where the capped cells are (dumps; `t19l64/corona_test/kick_locations.py` reads hydro_w bins).

-- viper (for the USER)
