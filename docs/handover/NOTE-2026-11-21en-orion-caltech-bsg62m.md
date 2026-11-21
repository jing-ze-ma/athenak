# NOTE 2026-11-21en orion -> caltech (USER): sync 2 new BSG tables + M_bsg62m (BSG production candidate) at the HEAD

1. **Sync your FRESH copy of fresh_prod_1120** from 18006c14b. There are two new BSG tables (plain text, not gzipped):
   - bsg/rosseland_mesakap_gs98_x0.70_z0.008.txt (808275 B, md5 7fa97c20db4cbe151279fc4f219605c1)
   - bsg/planck_ext2_gs98_x0.70_z0.008.txt (792760 B, md5 92fd7b56cb2892961afc61661fc42610)
   The deck below names them through FRESH_DIR/bsg/..., and armg.sh checks that every input exists.
2. **HEAD**, ahead of G_bsg62 / GR_bsg62 and after any arm already running: this is the BSG production candidate (MESA Rosseland matching Linhao's profile kappa to 1.4 %, mu 0.62, WB, THINF, maxit 800, conv_max 1e-5, tfloor 1000 K). Run on B_v203, 1 node, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys:
   - `RT=1 armg.sh M_bsg62m $B_v203 decks/bsg62m_wb_smoke.athinput 120` (deck in 18006c14b)
Please commit res / hst / log tails as soon as it finishes. Production still waits for the explicit go NOTE.

-- orion (for the USER)
