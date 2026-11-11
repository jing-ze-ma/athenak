# TASK 2026-11-11 viper -> ORION: vgdh-sw -- a wider thin-switch smoothstep (implicit_thin_h0 / _h1): cs order + no-regression

## Why (viper CPU, pulse_cs C R3 32-64-128, CAND + keys; /viper/ptmp2/jinma/rc18_1110/VGDH.md "g/w convergence")
The x1-face blend weight w = S(g) converges at p 0.91 (g itself 1.44) and multiplies an O(1), resolution-independent
branch mismatch (F_hr - F_c)/cE = -0.12 .. -0.06 in the band at the front (same at n 32 and 64). A wider smoothstep lowers
dw/dg and lifts the order with EXISTING keys (no code):
| h0..h1 | E p | F1 p |
|---|---|---|
| 0.3..0.6 (CAND) | 0.91 | 1.11 |
| 0.2..0.8 (sw2) | 1.64 | 1.70 |
| 0.0..1.0 (sw3) | 1.85 | 1.81 |
| 0.1..1.5 (sw4) | 1.90 | 1.80 |
w = 0.5 constant: 2.58; w = 1: 2.84 but breaks the lat twins (your 71d0c0bc0). sw3 / sw4 keep w -> 0 in thick cells
(g -> 0) but sw4 never reaches w = 1 at g = 1 (S = 0.77): the beam decks decide.

## Arms (keys ADDED to CAND; binary 92d4621fd for cs, fork/vgdlay-1110 453f879bd for sp)
- CAND-sw2: `implicit_thin_h0 = 0.2`, `implicit_thin_h1 = 0.8`
- CAND-sw3: `implicit_thin_h0 = 0.0`, `implicit_thin_h1 = 1.0`
- CAND-sw4: `implicit_thin_h0 = 0.1`, `implicit_thin_h1 = 1.5`

## Runs, in this order
1. No-regression on latmid_cs / latatm_cs (T and C, Amode + Emean/F1mean) and lat_th_atm / lat_th_mid: sw3, sw4 (then sw2).
2. Cartesian beam suite (docs/handover/beams-1110, k1_redist 0) on top of the Cartesian cand: sw2, sw3, sw4.
3. pulse_cs C R3 32..256 + T 8..64: sw3, sw4.  pulse_u SP_FULL (one fixed layout): sw3, sw4 (viper runs 64..256 too).
Build line as in NOTE-2026-11-08-orion-viper-consfix-ledger.md. NOTE back (NOTE-2026-11-11-orion-vgdh-sw.md) with the
table: arm x (pulse_cs C/T, pulse_u, 4 lat twins, cb/ws66/shd/hcyl).
