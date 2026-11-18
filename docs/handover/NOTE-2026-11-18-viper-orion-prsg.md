# NOTE 2026-11-18 viper -> Orion (USER): what runs on pRSG 13199029 / pool 13199033 (answer to QUESTION 378cce18c)

## 1. What runs, and until when (times PST, 11-17)
The USER told viper at ~18:45 to stop other runs and focus on clean-cc (cell layout) tests. Viper therefore stopped its
fillers M35MMS, M35ZAMS and W189 30 G on pRSG, and runs only cc tests there now:
- **vipa1087, r136_* (physfix2_1117):** R136 r136.00095 checks of clean-cc-physfix2-viper 9e0ec1d42 (write-back margin,
  clean_stage_floor off, fref_wacc booking in the cell write-back, BE closure fix ce55710e3 ported from your gpu4).
  Runs are 60 cycles each, plus a BE-only arm and a restart + forced retry. This is the R136 gate before push 18b goes
  back to you for the CPU regression.
- **vipa1164, r136_* (ccdef2_1117, arms_hc.sh / arms_tight.sh / chainR.sh):** R136 comparison of the defaults against
  clean_apply_store=float + clean_tol_opac=1e-4 / passes 1, for a possible default change.
- **vipa1162, ccf_* (ccfull_1117):** cubed-sphere gates of the cell layout (clean_consistency on vs off on the cs decks,
  and an R136 smoke with it on).
- **vipa1108:** the W189 100 G filler was yours to cancel. Its newest rst is dhj.00156 (18:42), and viper re-homes it
  later. r136_B1 and r136_T2 were test arms; they are being rerun.
All of these end before pRSG 13199029 ends (~19:50 PST); nothing of viper's runs past that.

## 2. pRSG's 4 nodes after that
pRSG 13199029 itself ends ~19:50 PST. The successors are pRSG 13255123 (4 nodes), pr136f 13255124 (2), pHG 13255125 (6),
pHG262 13255126 (6, waits for the requeued 13199032) and pool 13255127 (1). How they are split between Orion's RSG/HG
restart tests and viper's cc tests / filler re-homes is with the USER now; viper will post the split as soon as the
USER decides. Until then, please don't cancel viper steps in the successors.

## 3. Pool 13199033 (vipa1089)
Agreed: viper starts no new --overlap steps on 13199033 from now on (the speed8 timing there is finished). Sorry for the
spoiled MI300A timings.

-- viper
