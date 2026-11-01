# TASK 2026-11-01 viper -> Caltech (USER): WASP-121b A/B of the cube-vertex EMF key <mesh>/cs_seam_fc_nolimit, branch fork/csvertex-test-1101 @ 9f9a3d508

Reply to NOTE-...-hst25-wbtan (8d622b6bc).
- **Your item 2:** no viper run uses roche_wb_tan; viper runs no production. Thanks for the fix.
- **Your item 1:** thanks for the hst check.

## What the key does (default off = bitwise; not for fin51-r-fix yet)
- **The problem:** the resistive EMF (Ohmic and ambipolar) converges only at ~1st order in its maximum at the 8 cube vertices.
- **Cause:** the monotonicity clamp in the face-centred seam resample (bvals_fc.cpp seamval). On smooth data it fires and puts an O(h^2) error into one ghost face beyond each vertex. The two-pass curl turns that into an O(h) spike in J.
  - With exact inputs the vertex stencil is 2nd order, so the frame transforms and c/s factors are correct.
- **With the key on:**
  - the vertex x1e max converges at order 1.94-2.00 (nx2 16..256, Ohmic and AD alike), e.g. Ohmic 3.1e-3 -> 1.3e-5 vs key off 1.2e-2 -> 1.4e-3 (nx2 16 -> 128);
  - L1 is unchanged and div B stays at round-off.
- **Scope:** the key drops the clamp in the FC seam resample for ALL cs MHD (the field halo), not only resistive runs. The cell-centred resample keeps its clamp.
- **Viper gates:**
  - Key off bitwise vs b935e1452 on CPU: cs11 Ohmic 7/7 restart files; strong-field blast 4/4.
  - GPU: W121 s10/b100 restarts identical.
  - Key on, W121 s10, 169 cycles: finite; dt equal; E_B within 2.6e-8.
  - Discontinuity sweep: vertex and seam-midpoint blasts, beta 0.1 and 0.01, key off vs on. Key on runs wherever key off runs; |B| and p extrema, event counts, symmetry and div B agree to <= 0.04 %.
- **Branch:** 9f9a3d508 = b935e1452 + ambipolar-1101 (the merge afa7cca52, default off) + the key + cs_test diagnostics.

## Arms (r1m WASP-121b cs MHD deck with its production keys; same binary and same restart for both; 1-2 rot)
- **A:** control, no new key.
- **B:** A + `mesh/cs_seam_fc_nolimit=true` (adding it on the command line works on a restart; the log prints a NOTE).
- **Optional C:** A with a different rank count or node layout, as a chaos-spread reference. Otherwise use the panel-to-panel scatter.

Outputs: hst at the usual cadence; event log (file_type = log, dcycle 100); a bin and rst dump at 0.5, 1 and 2 rot.

## Diagnostics (A vs B)
1. Finite, no FATAL. Cost per cycle must match within 1 % (the key adds no work).
2. Summed event counts: dfloor, efloor, tfloor, vceil, c2p fail, maxit_c2p, fofc.
3. hst: mean and min dt, E_B, KE, total-E drift.
4. From the dumps at p > 1e-6 bar, in three regions (within 2 angular cells of a cube vertex; within 2 cells of a seam, excluding vertices; panel interior):
   - rms |B|;
   - rms |J|;
   - grid-scale power = rms(f - 3x3 angular box mean)/rms f, for B_r, B_t and v_r.

   Give the vertex/interior and seam/interior ratios for A and B.

## PASS
- B is finite wherever A is.
- B event counts <= 1.2x A (or within the A/C spread).
- B mean dt within 2 % of A.
- |E_B(B) - E_B(A)| and |KE(B) - KE(A)| within the spread.
- Vertex/interior and seam/interior grid-scale and |J| ratios: B <= 1.1x A. B is expected to be lower at the vertices.

The USER decides the default after this test.
