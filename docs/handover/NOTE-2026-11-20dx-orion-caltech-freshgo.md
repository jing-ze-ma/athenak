# NOTE 2026-11-20dx orion -> caltech (USER): FRESHGO smoke + block-shape arms (B_v203), HEAD of the queue

Files are in `docs/handover/fresh_prod_1120/freshgo/`: armg.sh (the runner), resf.py (the bitwise check), decks/*.athinput (6 fresh + 6 smoke decks) and **ARMS_smoke.txt** (the arm list with exact commands).
- **Build:** B_v203 (caltech-runs/v203 e688cb7bc).
- **Inputs:** set FRESH to a **gunzipped** copy of `docs/handover/fresh_prod_1120/` (all *.gz gunzipped in place). armg.sh substitutes FRESH_DIR in the deck and checks that every input file exists.
- **Runner:** `HOLDER=<job> NODE=<host> OUT=<root>/freshgo FRESH=<dir> MEMWRAP=<your memwrap.sh> [RT=0] armg.sh TAG $B_v203 decks/<deck> NC [keys]`, 1 node per arm (NR=4).
  - Several single-node arms inside a multi-node holder are fine: armg.sh uses `-w $NODE -N 1`.
  - **Always pass TAG** (armg.sh does `rm -rf $OUT/$TAG`).
- **Smoke arms** (a = 50 cycles straight; b = 25, rst, restart to 50; resf.py compare), in this order: S_t9l62, S_bsg, S_t19l64, S_m35mms, S_t9l62_corona, S_t19l64_corona.
- **Timing arms** (RT=0, 8 cycles): B4_t9l62, B16_t9l62, B8_bsg, B16_bsg, B8_m35mms, B16_m35mms.
- **Priority:** right after R_LF and the v20.3 gate. They **preempt CONVERR / LF seeds / EW** (logged for later restart). Run as many concurrently as nodes allow.
- **Please commit per arm:** res.txt, info.txt, the last 200 lines of a/run.log (and b/run*.log), a/*.hst and the mem peak, to `docs/handover/caltech_results/<TAG>/`. Bins can be deleted.
- Fresh **production** is still prepared only. It launches only on orion's explicit go NOTE.

-- orion (for the USER)
