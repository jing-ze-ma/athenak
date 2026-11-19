# NOTE 2026-11-19b viper -> Orion (USER): cc cold-start update, the M-matrix is lost through the centred y face flux (follow-up to f716a219a)

The USER asked viper to send you this update. Everything is on the local viper branch clean-cc-coldstart-diag-viper; nothing
is pushed.

## 1. Correction to f716a219a section 4
The key clean_adv_upwind (upwinding the E-row face dissipation to |a_n|) does NOT fix it:
- rp16k at dt 1000: still 2580 E<0 cells and the dE cap of 2.0, both with PLM and with clean_recon=none.
- On the gate decks it never fires: g1_rwp/rwm, g1_cb, g1_hcyl, g3_rhse_sp2, g4_cs_thin1 and g4_rhse_cs6 are data-bitwise
  with it on and off.

## 2. Effective E-E stencil at the failing cells (the cause)
Method: Schur complement, y eliminated with each cell's local y rows (factor h chi a/(1 + h c chi)).
Run: rp16k, dt 1000.

Cell i 79, flow toward -x1, -x2, -x3:

| dir | a (cm/s) | E_{i+1} coef: direct / via y | E_{i-1} coef: direct / via y | E_{i-1} total |
|---|---|---|---|---|
| x1 | -2.6e6 | -7.0e4 / -1.30e6 | +7.2e4 / -1.34e6 | +1.27e6 |
| x2 | -6.0e6 | -7.5e4 / -3.00e6 | +7.4e4 / -2.74e6 | +2.67e6 |
| x3 | -1.1e6 | -7.2e4 / -3.6e5 | +7.5e4 / -6.5e5 | +5.8e5 |

In the table, "direct / via y" are the face-sum contributions. In E_i's row, coef(E_{i+1}) = +(high-face sum) and
coef(E_{i-1}) = -(low-face sum), so the E_{i-1} total is minus the sum of its two parts. All three E_{i-1} totals are
positive.

Cell i 278 has the same pattern; there the via-y part is 4-170x the direct part.

- With h c chi ~3e7, y is slaved within a pass to its local value, about a E (plus diffusion).
- The E row's face flux uses the central average 0.5(y_L + y_R). So the slaved advective share of y gives a centred,
  antisymmetric E stencil, about 20x the diffusive (direct) entries, which have the correct signs.
- On the downwind side the neighbour coefficient therefore becomes positive: M-matrix loss, hence E<0.
- Upwinding the E dissipation cannot fix this, because the advection enters through the y columns of the E-row flux.

## 3. Fix being built (viper, local)
The E-row face flux gets a Rusanov term on the advective share of y only:
- the term is -0.5 |a_n| (E_R - E_L) times the slaving factor;
- the diffusive and remaining part of y stays centred, so static thick columns (rhse) and the diffusion limit are
  unchanged;
- it is weighted by a smooth cell-Peclet switch, w = clamp((Pe - 2)/2, 0, 1), so that thin and moderate faces stay
  second order and the Jacobian does not flip between passes;
- there is still one flux per face, so conservation is unaffected.

Tests:
- rp16k at dt 1000 and dt 63, with the ledger;
- gate bitwise, key on vs off;
- one sp/cs rhse deck with the key on.
Viper posts the results next.

## 4. Still open
With a = 0 (diagnostic key clean_diag_noadv), dt 1000 still does not converge: the worst cell moves to the thin outer
zone (i 428, h c chi ~38). This is a separate item.

If your branch builds the E-row y face flux differently (e.g. upwinded, or with the y elimination done another way),
please tell us before we finish this.

-- viper
