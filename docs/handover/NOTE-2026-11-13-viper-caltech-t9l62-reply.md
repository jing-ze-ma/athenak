# NOTE 11-13 viper->Caltech: reply to e393b895e (T9L6.2 seg 1)

- Restart: viper GPU fresh-vs-fresh IS bitwise (MI300A, same deck, 512x64x128). Restarts are bitwise EXCEPT when the
  rst lands right after a stage fallback + BE redo: the next (auto_be skip) BE step differs at round-off. We are bisecting
  the missing state now (not deck-specific; it probably affects any hesdirk2 + auto_be run). Fix + gate will be pushed on
  a fork branch. Until then continuing from t9l62.00010 is fine physically, just not bitwise.
- Density floor at cycle 0: expected-ish. The IC (Jiang init_outburst) is the star; above its outer edge the column is
  capped at the floor, and the deck's x1max is beyond the IC edge. Please report the radius where the floor starts at
  cycle 0 and whether it fills in (Mdot_top) by 111 h.
- TOI-1518b 20 G runs on Raven (UPDATE2 504da16cb); nothing dhj is needed at Caltech, ignore the bundle.
- Continue T9L6.2 from t9l62.00010 in the next 2-GPU slot when Slurm is back; send L(r_out,t)/Mdot to 111 h.
