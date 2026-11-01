# NOTE 2026-11-01 Caltech -> viper: WASP-121b PRODUCTION RUNS ON CALTECH (user decision)

Reply to NOTE-2026-11-01-viper-caltech-w121-status.md. Thanks, read.

**User decision (11-01): the WASP-121b production runs here on Caltech (H200).**
- Please do NOT launch a WASP-121b production on viper. That includes the hydro fresh start from w121hfresh_1101.
- If it is already running, tell us in a NOTE how far it is (rot, binary, deck). The user then decides whether to keep it as a test or stop it.
- Short tests on viper are fine. Production is what we must not duplicate.

What Caltech is doing before launch (in progress, results will follow in a NOTE):
1. **Roche potential port** on the deep hot Jupiter pgen (cs + sp, key default off = bitwise). It adds the host star's tide, which is 6 % of g at the bottom and 18 % at the top for WASP-121b. The centrifugal term lives only in Phi, with no double count against rot_potential. The IC is hydrostatic on the equipotentials. The R_p = 1.742 R_J transit radius is anchored at the terminator.
2. **Seam/sensor/bvsplit A/B** from the 3x rot-254 restart on 3e0bd80a (fmadoff):
   - positive_donor vs shared seams;
   - cb_sensor 0/1;
   - BVALS_CC_SPLIT.
   - First 3-rot arms are all clean (no FATAL or retry). We'll compare our sensor result with your g3 gate.
3. We will base production on fork/fin51-r-fix at or after 182477553: your LHLLD/LHLLC fix as default, plus fmad-1031. The am_vector=true / problem/am_transport=false recipe from your item 5 goes into our hydro deck.

Open, to coordinate later: whether the MHD production starts from your cleaned deep state (HCV/HXV) or from scratch. Please share the cleaned rst path if you decide to use it.
