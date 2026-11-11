# NOTE 2026-11-11 ORION -> viper: vgdh / vgdh-w1 sp follow-up -- pulse_u full refinement on 453f879bd (one block per level):
# col / lf0 / lf2 / wsx0 barely move CAND (finer pair 0.45-0.62); w1 diverges; only the wider smoothstep helps (sw3/sw4)

The sp rows promised in 978eeecf6 (vgdh), 99e55db19 (vgdh-geo) and 9f868b843 (vgdh-w1). fork/vgdlay-1110 **453f879bd**
(the vgd_hr ghost fix), my st17 SP_FULL: nx2 = nx3 = max(4, nx1/16) = 4 8 16 32 for nx1 64 128 256 512, ONE angular
MeshBlock per level (SP_MBN=1: layout fixed in physical angle, as addendum 1 asked; 1 rank per run). Rotation off, CAND +
arm keys. All levels rc 0, no FATAL. (My earlier radial-only / 4x4-block pulse_u numbers are void, as you said.)

| arm | E p (64-128-256, 128-256-512) | F1 p |
|---|---|---|
| CAND | 1.14 0.45 | 0.87 0.56 |
| CAND-col | 1.10 0.62 | 0.85 0.73 |
| CAND-lf0 | 1.09 0.48 | 0.81 0.59 |
| CAND-lf2 | 1.12 0.48 | 0.91 0.58 |
| CAND-wsx0 (radial faces central) | 1.24 0.61 | 1.02 0.80 |
| CAND-w1 | 1.67 **-1.23** | 0.96 **-0.99** |
| CAND-wc5 | 0.73 1.44 | 0.39 1.90 |
| CAND-sw3 (0.0 / 1.0) | **1.31 1.53** | 1.05 **1.73** |
| CAND-sw4 (0.1 / 1.5) | **1.31 1.41** | 1.16 **1.92** |
CAND-geo3 / geo4 and CAND-sw5 still running (follow-up). Cost at C_256 (1 rank, 285 cycles): ~750-880 s per run for the
geo-2 arms.

## Reading
Unlike cs (where wsx0 = central radial faces gives 1.91 1.92), on sp even wsx0 stays ~0.6 at the finer pair: the sp
residual is not (only) the radial half-range faces. The wider smoothstep is the one change that lifts sp (E 0.45 -> 1.4-1.5,
F1 0.56 -> 1.7-1.9), consistent with the cs finding that the radial variation of w = S(g) is the O(dx); sp would need the
next level (nx1 1024, angular 64) to tell pre-asymptotic from a further cap -- not run (1 rank per run, ~hours).
