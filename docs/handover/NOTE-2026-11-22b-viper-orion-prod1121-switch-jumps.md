# NOTE 2026-11-22b viper -> orion (USER): prod1121 HG and RSG jump at the restart into the clean staggered scheme

We read these from the page reductions of your prod1121 segments against the last old-scheme dumps. They are numbers only,
not diagnosed. The pages show them with a red switch line and a "recent" history figure: HG
https://claude.ai/artifact/Twfp2mnaegk4asaJXPHdfg (v10), RSG https://claude.ai/artifact/5zYUnhQSYYjt39vT1XMsNp (v14).

**HG** (restart hegiant.00201, 760.42 d; Raven `stag2/prod1121/hg`, dumps 02282-02286 to 762.0 d):
- **Cells at the 2400 K T floor:** about 390 in old dumps 02280-02281. Every one of your dumps has 13.7-17.2 M (about
  20-25 % of the mesh).
- **L_bot/L_in:** about 1.05 -> about 1.5.
- **Median R_ph:** 0.95 -> 1.04 of the IC radius at the switch.
- **Mdot_top:** up by about a factor of 100.
- **Health:** FATAL 0, NaN 0, non-converged 0; 3.31 s/cycle.

**RSG** (restart rsg.01707, 8361.09 d; viper `orion_stag2/prod1121/rsg`, dumps 08362-08391):
- **L_mid/L_in:** about 0.15 -> about 0.52.
- **L_top and L_bot:** spike briefly, then settle near 1.
- **Median T_eff:** about 3550 K -> about 3455 K.
- **Health:** FATAL 0, NaN 0, 0 non-finite cells; 2.24 s/cycle.

**R136:** the page is still being reduced; I will add it if it shows the same.

These are jumps at the restart itself, not slow drifts. The HG floor fraction stands out most. Please check whether the old-
to-clean restart carries the radiation state over consistently: the m1/vet state, the staggered layout fill, and the T
floor, or keys such as floorcap and the floor-gas switch that are on in prod1121 but were not in the old decks. The full
old key sets are in `prod_1107/hg/run.log.13199032.0` and `prod_1107/rsg/run.log.13359763.0`.

-- viper (for the USER)
