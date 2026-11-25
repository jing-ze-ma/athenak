# NOTE 2026-11-24zzzzc4 orion -> Caltech: t9l62, which cells sit at the velocity ceiling (USER: "ask caltech to check")

Read-only analysis of files that EXIST; no simulation, nothing written into a production directory, no key changed. CPU on an idle or shared node; low priority against production, t9l62's start and the c4 items. vceil = 3.0e8 cm/s in the t9l62 deck.
1. Dumps / rsts: any t9l62 files you hold with 4.6e5 < t < 5.7e5 s (outburst), 3 to 6 times spread over it, plus one quiet (c2553 stop rst or the latest quiet) and one late (00037 or newer, old scheme: mark it). Name each file, its run, binary and t. If a needed time has no bin but has a rst, read the rst. If none is readable without a run, say so and stop.
2. Per file, over the real cells (no ghosts): number, volume, mass and kinetic energy of cells with |v| >= 0.5, 0.8, 0.9, 0.99 and 0.999 vceil; |v| maximum with its cell (r, theta, phi, rho, T, v components); mass and KE as a fraction of (i) the star above r = 2e12 cm and (ii) the mass and KE of all cells with rho > 1e-10. Also the radius range, rho range, tau_cell, and c_g / v distribution of the >= 0.9 set. A histogram of |v| in the rho < 1e-11 gas (log bins 1e7 to 3e8).
3. Are they clipped? A cell AT vceil to printed digits (|v| within 1e-6 of 3e8) is a clip; a smooth tail below it is not. Give both counts. State separately what the dumps cannot show: the ceiling acts in c2p and its energy goes to heat (vceil KE -> heat) inside a step.
4. From existing logs of t9l62 runs (run.log of S5G, T2, T9 arms, c4 start): the per-cycle `skip(dfloor/vceil)` and `vceil` / FloorGasCap counters (floor_audit lines), how many cells and which fraction are velocity-ceiling events as opposed to dfloor; say if the counter mixes the two and can't be split (then don't split it).
5. Return: a table per file, the histogram as text, the counter table, and file list. No interpretation needed beyond what the numbers show; I score it.

-- orion (coordinator)
