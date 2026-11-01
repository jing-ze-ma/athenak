# NOTE 2026-11-01 viper -> Caltech (USER): skip G4 / rkghalo (part B of 6e5e226a6); viper measures the new version

- **v1 (0594b145d) is bitwise but useless:** it falls back to the full exchange in ~99 % of super-stages and is 2.4 % slower. Drop it.
- **v2 (fork/rkghalo-test-1101 @ 3503c1661):** exchanges only the changed slots, which is {IDN, IEN} on every stage of the s10 state.
- Viper is gating v2 (bitwise, cost, profile on apudev). It only matters for non-lean runs, so no Caltech GPUs are needed.
- **Thanks for the csvertex A/B.** The default (fresh starts only) is being merged into fin51-r-fix now; we will NOTE the hash.
- **max_eta:** waiting for your final R14 at +1 / +1.5 rot.
