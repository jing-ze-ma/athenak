# STITCH.md (10-09): ext2 Planck + orion LTE low-density shape (AG Car, X 0.36, Z 0.02, GS98)

**Table:** `planck_ext2orion_gs98_x0.36_z0.02.txt` (md5 493f9d01). It replaces `planck_ext2_gs98_x0.36_z0.02.txt` (dc482452, branch agcar-opac-1009).
**Rosseland:** not changed. Use `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5) from agcar-opac-1009 as it is. It is not copied here.

**Format.** The format is the ext2 format: same 91370 lines, comment lines at the same positions, grid line `# 217 421 2.6 0.025 -21 0.05`, T-slowest, `%.5f`. Two header lines were reworded, with no line added:
- line 3 now describes the stitch;
- line 9 (VALID BOX) has an appended clause.

The reader is `HsReadOpacityTable` (src/pgen/he_star_m1.cpp l. 218), as `he_planck_table`, which is what TABLES_EXT.md names. `read_table()` in `stitch_ext2_orion.py` mirrors it: comment lines are skipped, the first comment that parses as `int int real real real real` is the grid, and the code asserts value count = nT*nD. The new table parses: 217*421 values.

**Rebuild:** `python3 stitch_ext2_orion.py <agcar-opac-1009 dir> <agcar-opac-orion-1009/tables> <outdir>`, then `check_stitch.py <outdir> <ext2 planck>`. Outputs: `checks.txt`, `floors.txt`, the two png files.

## Recipe
1. **Floor rho_f(T)** = the lowest ext2 node that holds real source data. It is reproduced from `build_ext.py fergR` (`floors.txt`):
   - log T <= 4.0: the Ferguson 2005 box, lD - 3 lT + 18 >= -8. This is the script's own `RRg >= rs[0]` mask, inside which it overwrites the rows with `bilin` of the data. rho_f runs from -15.65 (3.45) to -14.00 (4.0).
   - log T >= 4.2: TOPS. Per TOPS row the clamp floor is found with the script's scan (min with -14). A node is real only if both TOPS rows bracketing it in log T are real there. rho_f = -15.00 at 4.20-4.30, from the rows 4.1615 (-15.2), 4.2407 (-15.0) and 4.3657 (-15.0).
   - Blend 4.025-4.175: max of the two floors. This is the Ferguson box, -13.90 to -13.45. rho_f therefore drops by 1.55 dex between 4.175 and 4.20, which is ext2's own source structure.
   - Check: below rho_f, ext2 is the power law. Its max second difference in log rho is 1e-5 (rounding) on the Ferguson rows and <= 0.009 on the blend and TOPS rows, where the two bracketing TOPS rows have different floors. The kink sits at rho_f.
2. For rho < rho_f: log kP_new = log kP_ext2(T, rho_f) + [log kP_orion(T, rho) - log kP_orion(T, rho_f)]. Cells at rho >= rho_f keep the original text.
3. **T edges.** The correction delta = new - ext2 is multiplied by:
   - 1/3 at the orion edge nodes log T 3.450 and 4.300;
   - 2/3 at 3.475 and 4.275;
   - 1 in between;
   - 0 outside.

   This is a linear ramp from 0 at the first node outside (3.425, 4.325) to 1 at the third node inside. Orion starts at log rho -21, which is ext2's lowest node, so no hold below -21 was needed. ext2 has nothing below -21.
4. Rosseland: unchanged (see above).

## Checks (`checks.txt`)
- **(a) Cells changed:** 4470, all at log T 3.45-4.30 and rho < rho_f. Max |delta| per region:

  | region | cells | max abs delta (dex) | mean delta (dex) |
  |---|---|---|---|
  | taper 3.450-3.475 | 216 | 0.131 | +0.013 |
  | Ferguson 3.500-4.000 | 2630 | 1.616 | -0.019 |
  | blend 4.025-4.175 | 1024 | 0.997 | +0.296 |
  | TOPS 4.200-4.250 | 360 | 0.495 | +0.198 |
  | taper 4.275-4.300 | 240 | 0.410 | -0.065 |

  At log rho -21, delta runs from -1.62 (log T ~3.7, where H- falls away; orion is lower) to +1.11 (log T ~3.9, lines; orion is higher).
- **(b) Continuity across rho_f** (|kP(rho_f) - kP(one node below)|): max 0.018 in ext2 and 0.039 in the new table (log T 3.850). The stitch is continuous in value by construction. The 0.039 is orion's own slope over one 0.05-dex node.
- **(b) T-direction** (max |row-to-row jump| at log rho below 13.45, log T 3.40-4.35):
  - overall: ext2 0.570 (3.950-3.975, -21), new 0.683 (3.650-3.675, -21);
  - low edge 3.40-3.50: ext2 0.202, new 0.227;
  - high edge 4.25-4.35: ext2 0.106, new 0.336.

  The largest high-side jump, 0.469, falls between 4.225 and 4.250 (full weight, not in the taper), against ext2's 0.175 there. It is orion's own T structure at -21. The untapered correction at the edge rows is max 0.098 at 3.450 and 0.534 at 4.300; the taper spreads the latter.
- **(c)** Unchanged cells: 0 text differences among the 86887 cells with delta = 0. Data lines that differ = 4470 = the changed cells. Written vs computed: 3.3e-6 (`%.5f` rounding).
- **(d) Plots:**
  - `stitch_delta_map.png`: delta over log T x log rho [g/cm^3], with rho_f(T);
  - `stitch_cuts.png`: kP(rho) at 4000, 5000, 6000, 8000, 10000 and 15000 K (nearest nodes 3.600, 3.700, 3.775, 3.900, 4.000, 4.175), showing ext2, orion and the stitched table, with rho_f marked.

CPU only (login node). No AthenaK source change, no restart used.
