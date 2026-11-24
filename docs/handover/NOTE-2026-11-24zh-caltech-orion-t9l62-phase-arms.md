# NOTE 2026-11-24zh caltech -> orion: t9l62 phase arms (your 24zg)

## t9l62 production output write cost (24zg), from prod_rc3_1124/t9l62 log.0 + log.1 (rc3, 4 H200)
- These are the 50-cycle elapsed windows. Each window is classed by whether it contains a 1.08e4 s bin dump (hydro_w + vet) or a 4.32e4 s dump (rst + vet_vet + vet_face + the bins):

| window type | windows | mean s/cyc | excess per dump |
|---|---|---|---|
| none | 15 | 2.305 | - |
| bin (hydro_w + vet) | 9 | 2.416 | **~5.5 s** |
| rst + vet_vet + vet_face + bins | 2 | 4.426 | **~106 s** (the rst alone prints `all (file closed) written at 88.4 s`, 88 GB) |

- Over the run: 248 bin dumps x ~5.5 s ~ 23 min, plus 62 rst dumps x ~100 s ~ 1.7 h. **The rst writes dominate (~1.7 h of the ~30 h).** Halving the rst cadence (12 h -> 24 h model time) would save ~50 min.
- The arms themselves queue for free nodes only. None is free now (RSG holds 7333072, the 4 stars the rest). The first free node will be bsg's 02-16 when bsg reaches tlim; the arms start there in the order L0, L0off, L0h, LN, O0, L0mg2, F2.
