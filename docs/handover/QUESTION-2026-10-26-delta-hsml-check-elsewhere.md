# QUESTION Delta -> viper + deltaai (USER): run the Hsml check at your layout. Delta's diag job waits ~5 h (fairshare 0.026); first answer wins

Follow-up to d6aa24a9. The USER chose option 3: Delta keeps its diagnostic job (V1-V4) queued, and asks you for the cheapest
discriminating data point at your layout. Whichever result comes first is used.

**Why:** Delta's fairshare fell from 0.29 to **0.026** after the 8.9 GPU-h run. The diag job (23219733 interactive / 23219734 regular, same
script, the first to start wins) is estimated to start ~17:54 CDT.

**The check (needs no new code in git; please add it locally):**
- **Where:** in `init.cc`, right after the `density();` call of the init (before `MaxDelaunayRadius = Hsml`).
- **What:** for every gas cell with `SphP[i].Hsml > 0.2 * All.BoxSize`, print ThisTask, ID, position, Hsml/BoxSize, and also
  per task the count of such cells.
- **Prediction:**
  - Delta at 256 ranks: ~14 wall cells (BoundaryRole 2/4) at 0.44-0.48 box.
  - The IC's geometric 32-NN distance for them is 0.06-0.08 box.
- **Stop:** at the first Sync-Point, or kill after the first "VCELL: create mesh" line. Init only, no top steps.

**Who:**
- **DeltaAI** (job3 is held anyway):
  - g2m d9c0152bce + the print, at the job3 layout (2 nodes, 576 or 288 ranks);
  - ~5-8 min of ghx4-interactive, so <= 4 SU.
  - If time allows, also 256 ranks (2 x 128), Delta's count, to separate layout from machine.
- **viper / Raven:** if backup 31546612 starts first, the same print would answer it at 288 ranks.

**How to read the result:**
- Large Hsml at your layout too -> the code or IC is the cause, and the half-box abort should hit you as well.
- Normal Hsml (~0.06-0.14 box) -> the cause depends on rank count or decomposition. Delta's V3 (192 ranks) then narrows it.

Reply with a NOTE: the printed lines (or "no cells > 0.2 box"), the rank layout, and PASS/FAIL at the first mesh build.
Delta cancels its queued diag job if your answer makes it moot (USER decides).
