# NOTE 2026-10-23 (viper): implicit_precond_fallback merged into fin51-r-fix (c1681dbe3)

**fork/fin51-r-fix is now c1681dbe3.** It was fast-forwarded from d01ce628c.

The new commit adds one opt-in key, off by default:
```
<rad_m1>
implicit_precond_fallback = none    # none (default) | rbgs | line
```

## What the key does
A BiCGStab inner solve that stops at `implicit_lin_maxit` with its stopping test unmet is redone once with the fallback preconditioner. The redo:
- starts from the same initial guess;
- uses the same Eisenstat-Walker stopping test.

A breakdown does not trigger the redo.

When the key is set, the bicgstab summary line prints the number of redone solves and their inner-iteration total. Restrictions:
- rbgs as the fallback needs `krylov_fuse >= 1`;
- the key is refused on cs.

## Why we added it
On the IDORT vacuum-beam deck cbh (1024x256), `precond = mg` fails as follows:
- Almost every mg inner solve hits the 200-iteration cap, because coarse-grid corrections spread error into near-vacuum cells (E 1e-30..1e-20).
- The per-cell |dE|/E Picard test can never be met there, so the solves are reported NON-CONVERGED.
- That triggers dt_retry, then the BE redo, then FATAL in PROD.

The solution itself is unchanged by the switch: E with rbgs or line agrees with mg to 1e-21.

## Results
cbh, GPU, nlim 20:

| precond | NC | t reached | E_L1 | GPU s |
|---|---|---|---|---|
| mg | 11 | 0.027 | - | 426 |
| mg + fallback=rbgs | 0 | 0.080 | 0.580011 | 65 |
| rbgs | 0 | 0.080 | 0.580011 | 23 |

## Recommendation
- **Vacuum-beam / near-vacuum test decks** (cbh, cbg): set `implicit_precond = rbgs` directly. It is 2.8-7x cheaper than mg + fallback, because the first mg attempt runs all the way to the cap.
- **Production decks:** nothing changes.
  - The accretor D_B1v6 is bitwise identical to d01ce628c both with the key absent and with fallback=rbgs (0 redos; inner iterations mean 3.6, max 12).
  - You can set fallback=rbgs as a safety net. It costs nothing unless mg caps.
- The key must be present in the deck: a command-line-only key FATALs on a fresh `-i` run.

## Related side finding (no code change)
`implicit_lin_cnorm` 1e-10 (default) lets BiCGStab accept the initial guess when the per-step relative change in E is about 1e-10 or smaller.
- On a thin radiation wave at amplitude 1e-5 and nt=8N, that cost 4 % of the wave amplitude per period.
- Setting `lin_cnorm` 1e-13 fixes it.
- It only matters for tiny-amplitude, small-dt tests; production thin regions change far more per step.
- Check: a production log with "inner iterations mean" well below 1 would be suspect.

Details: viper /viper/ptmp2/jinma/blend_1022/PRECOND_FB_RESULTS.md, CBH_RESULTS.md, RESULTS_side.md
