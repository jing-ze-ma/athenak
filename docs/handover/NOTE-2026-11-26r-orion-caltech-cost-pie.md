# NOTE 2026-11-26r orion -> Caltech: module cost of the star runs (for a pie chart in the docs; the USER, 12:01 PST)

The USER wants a pie chart in the documentation like Fig. 20 of "ma2026", showing the fractional cost of each module (hydro, rays, rad-vet, others), representative of your star runs. Please profile them; Orion draws the figure.

**Arms** (one idle test node, 4 H200, 4 ranks; nothing on production nodes; NOT HG: nsys hangs there; not before 26m's suite is done): bsg, m35mms, m35zams, t9l62 (t9mlt only if cheap), each from a COPY of its newest rst, with its own production keys, binary athena_next3_c3f133b44 (md5 ad2f2007; the released code):
1. plain 40 cycles (no profiler): s/cyc, dt, passes / solves per cycle, stage rejections, memory peak per GPU;
2. then 25 cycles under your profc.sh tool (nsys phase only, NCU_ON=0), all 4 ranks, window without the first cycle: kname_r0..3.txt, gaps_r0..3.txt, stats. If nsys hangs on a star, stop it and use the Kokkos kernel timer (fenced; state the fence cost).

**Return** (text, no binary traces): per star, wall per cycle (plain), GPU busy and idle % of wall, and per-kernel ms per cycle (rank 0 and the mean of 4). Module map proposed (Orion will use fam.py; say if a kernel does not fit): hydro = hydro + gravity + EOS + history; rays = ray sweep, ring, GD face sums, closure from GD (families "closure"); rad-vet = Picard RHS / face fluxes, F solve, gas Newton / velocity / drag, E solve (operator + preconditioner), FCT, opacity, halo inside the radiation solve; others = the rest (boundary exchange outside, output, idle / host wait, unattributed). Give the sum of GPU-idle time as its own line, split MPI / sync / launch where gaps.txt allows. If the rays are refreshed not every cycle, give the per-cycle mean over the window AND the cycle with / without a ray rebuild.

One state per star (say early or late; a developed one, not the start-up); the wall time must be the one quoted in the docs table (26o item 2) for the same phase.

-- orion (coordinator)
