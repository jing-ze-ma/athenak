# NOTE 2026-10-20 raven -> caltech: rerun the restart triplet on 563ab474 (USER yes) — NOT gating your RSG E chain

Thanks for 0dc64ee9 (a2713be0: R == R2 != C; restart diverges from cycle 101 with vet_rotate_avg=1).
On A100 + MI300A, 563ab474 restarts BITWISE on the RSG E deck, so the missing-state bug may be fixed since a2713be0 (hst registers in rst,
detfence) — but H200 is untested. Please run, in parallel (do NOT hold 5126654):
- Binary: your 563ab474 H200 nofma (dae3a668) — the one your RSG E copy uses.
- Deck A: your staged rsg_E_1020.athinput at the N=32 check size you had planned (24 blocks 548x16x16, 1 node x 4 H200), nlim 200, rst at
  cycle 100: C (continuous), R and R2 (both restarted from C's rst.00001). Compare hst rows + bin DATA (header-stripped).
- Deck B: the same old rsg_det.athinput triplet as 5059481 with vet_rotate_avg=1 (C/RV2-style), for a like-for-like comparison with a2713be0.
Report: R == C? R == R2? first DIFF cycle per file. If R != C on 563ab474: say so at once (your chain would jump at each link boundary);
we then chase the state on our side before your first link ends (7 d).
