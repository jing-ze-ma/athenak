# NOTE 2026-09-29 (viper): opacity-Newton guard for <rad_m1> (m1-opn-guard)

For DeltaAI and Caltech. `implicit_opac_newton` has been DEFAULT ON since the m1-perf-0928 merge
(d26b7364; He box -58..-65 % ms/cycle). It has a defect. This merge adds a guard, which is ON
by default. Nothing needs to be done beyond pulling rt-integration.

## Defect

`M1OpnCell` (src/rad_m1/rad_m1_implicit.hpp) adds the Newton term of the flux opacity to the x1
rows as `coef += s c dt kE` and `rr -= s R_k`, with `s = q kt/B_k` and `q = -0.5 chat dt th G_f`.
The sign of `s` follows the face flux G_f and d kappa/dT (it flips across the Fe bump), and its
size has no bound. In fast, optically thin plume cells of the He presn wedge it reached
**-10035 on a diffusion diagonal of +802** (row diagonal TB = -9200; row-term decomposition from
he-presn-m1 f38bb18f, /viper/ptmp2/jinma/hepresn3d_0929/M1/diag_seam/D). The row is then not an
M-matrix: the solved E goes below zero, hesdirk2 stages fail, BE redos follow, floor clips are
applied silently, and the run eventually produces NaN. With Newton off (arm N) the solution up to
the failure is identical to 1e-8, so the term changes only the iteration path, not the fixed
point.

## Fix (src/rad_m1/rad_m1_implicit.cpp, the m1_impl_asm kernel; Cartesian and sp rows)

Each face's Newton contribution is recorded: the diagonal part, the neighbour-entry part and the
rr part. At the end of the row a face is dropped *entirely*, meaning all three parts in that row
(the face then stays Picard), when:

- it leaves the diagonal below `implicit_opac_newton_guard` (default **0.5**) times its value
  without that face; the faces are tested i+1/2 first, then i-1/2, so the diagonal can fall no
  lower than guard^2 of its value without either term; or
- (`implicit_opac_newton_guard_mode` bit 2, default **mode 2**) its neighbour entry ends up
  positive, which breaks the M-matrix.

Mode bit 1 (an rr test) exists but is NOT recommended: it drops almost every face (Picard 48
passes per solve on the wedge). `guard <= 0` gives the old unguarded rows exactly (bitwise).
Counter at the end of the run:
`<rad_m1> implicit_opac_newton_guard=...: Newton face terms dropped (all ranks, all passes)=N per solve=...`.

## Gates (all on the merged src; GPU apudev, 2 x MI300A, ROCm 7.2)

| gate | result |
|---|---|
| (a) He box bench (bench-2026-09-29-hebox), 50 cycles, guard=0 vs rt-integration tip 423bce0a | hst **bitwise** (job 12029873) |
| (a) dhj WASP-121b 1x bench, 50 cycles, branch vs tip | hst **bitwise** (job 12029874) |
| (b) He box, guard default | 0 faces dropped -> hst **bitwise** vs tip; Picard 2.12/solve (50 cycles) |
| (b) He box ms/cycle, 700 cycles, window 200-700, interleaved | guard=0: 40.05 / 39.63 / 39.70; default: 39.57 / 40.08 / 39.50; the speed-up is kept (Picard 2.009/solve in both) |
| (c) He presn wedge, ad3d rst 00016 -> t 19450, cfl 0.9, 2 ranks (he-presn-m1 + fix, arm F, job 12029875) | floor clips **0** (D: 463, N: 0); NON-CONVERGED **0** (D: 5); Picard mean **8.14** (N 10.25, D 19.45); faces dropped 704 (8.0 per solve); stage fallbacks 18 (N 16, D 22) |
| (c) variants | diagonal test only (mode 0): 7 clips, Picard 8.56; mode 0 with guard 0.9: 4 clips; mode 1 or 3: 5 / 0 clips at Picard 48; mode 2 with guard 0.9: 0 clips, Picard 8.23 |
| (d) tst rad_m1 (slab_cpu, opcheck_mpicpu, restart_mpicpu) | 3 passed; cpplint: no new violations |
| (e) hst at t 18992.49 (last row <= 19000), F vs N | max relative difference 5.3e-10 (solver tolerance 1e-8) |

Run dirs: /viper/ptmp2/jinma/opnguard_0929/{hebox,dhj,wedge/{F,G,G1,G2,G2t9,G3,G9}}.
The fix is also cherry-picked onto he-presn-m1 (local, not pushed).
