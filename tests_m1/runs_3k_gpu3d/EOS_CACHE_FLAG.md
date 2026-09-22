# The eos_cache "max |dq|/q" flag: 6.7e-2 on 1 rank (P1), 8.28e4 on 2 ranks (P2)

Diagnosed 2026-09-22 against HEAD 92d36db3. No code was changed.

## Verdict

The problem is in the diagnostic only. The solution is correct. The 2-rank value comes from
an ill-conditioned relative error, not from the cache, the MPI decomposition or the reduction.

## How the number is made (src/rad_m1/rad_m1_implicit.cpp, lines ~3739-3791)

- Once per step, for each interior cell (`m, ks..ke, js..je, is..ie`, local `m`, no halos),
  the check rebuilds the exchanged energy `q = emis*T^4 - dt*ch*kE*(de0+ep) + emis*4T^3*(r + ...)/b`
  twice: once with the cache values (`qc`) and once with the true table (`qt`).
- It stores `|qt-qc| / max(|qt|, 1e-300)` in the scratch slot `iw_(m,igm)`, takes a
  Kokkos Max, and then an `MPI_Allreduce(MPI_MAX)`. Both are correct.
  `ec_tmax = max over steps` is printed at the end (line 2445).
- The scratch slot `igm` doubles as the miss counter. That counter is summed before the check
  overwrites it (line 3716) and is zeroed at the start of every step (line 2600). So nothing
  the check writes reaches the state.
- The cache itself is per cell and per local `m` (`ec_(m,...,k,j,i)`), and is filled and read
  at the same cell. Nothing here depends on lid versus gid or on ghost cells.

The cause is the denominator. `q` is the net emission minus absorption. Where the gas and the
radiation are close to equilibrium (the optically thick interior), `q` comes from two large
terms that nearly cancel, and it passes through zero. There `|qt|` can be arbitrarily small,
while `|qt-qc|` stays at the cache round-off level (|de|/e = 2.5e-14 in both P1 and P2).
The largest ratio over about 5e5 cells x 600 steps is therefore set by the one cell that
lands closest to q = 0. Round-off (for example BiCGStab global sums over 2 ranks) moves that
cell. The range S 4.0e-3, P1 6.7e-2, PP 1.2e-3, P2 8.3e4 is the scatter of that extreme,
not a trend with rank count.

## Evidence that the solution is the same

1. GPU P1 vs P2 (bench/m1_gpu3d_0922/runs/P{1,2}/m1slab.{hydro,user}.hst). There are 91 rows
   at matching times up to t = 90.1. The largest relative difference is 3e-8 for tot-E,
   4e-9 for the KE columns and 1e-11 to 3e-7 for the user columns (F1top/F1bot/Etot).
   The 2-mom and 3-mom columns differ by up to 4.6e-2, but only on values that are near zero.
   The final dumps are at different cycles (648 vs 563, wall limit), so they cannot be
   compared directly.
2. CPU MPI reproducer: bench/m1_eosflag_0922/, snapshot of 92d36db3, gcc/14 + openmpi/5.0.
   It uses the same input with nx2 = nx3 = 16 (8^2 blocks) and 30 cycles.
   - The flag is 3.48e-4 on 1 rank and 3.30e-4 on 2 ranks (runs/c1, c2 run.log).
   - The final dumps (float32) agree to 1e-7 relative, which is one float32 ulp, in all
     variables and in all 4 pairings:
     cache 1 vs 2 ranks, no-cache 1 vs 2 ranks, and cache vs no-cache
     (`python3 bench/m1_eosflag_0922/cmp.py`).
     The large pointwise differences occur only in near-zero velocity components and F2/F3.

## Suggested fix (diagnostic only)

Normalise by the gross exchange rather than by the net:
`|qt-qc| / (emis*T^4 + dt*ch*kE*|de0+ep| + tiny)`.
Or report the absolute |dq| summed over cells relative to sum(|q|).
