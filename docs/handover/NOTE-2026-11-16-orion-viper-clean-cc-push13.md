# NOTE 2026-11-16 Orion -> viper: origin/clean-cc push 12 + 13 (62d71accc)

- Push 12 (0eeac5fdb): the FV-ray consistency-bracket positivity limit acts only in cells that went E < 0 in the
  positivity passes (it had limited steady starts at huge h). sp WB back to 6e-8..8e-8 (fv1/fv/lld).
- Push 13 (62d71accc): the cell-layout Rhie-Chow term uses its steady h-free coefficient (1-al) c/chi_f (was
  (1-al) h c^2/(1+h c chi_f), an O(dt) spatial-operator change): time order 2 restored in stiff radiation-heavy gas
  (P >= 0.03, tau >= 100: was 0.1-1). v = 0 decks change at truncation level; hcyl +0.4 %.
- Orion regression green on both (zero CLEAN WARNING, restarts bitwise). Neither touches the sp/cs GPU-divergent
  paths: TASK cc_gpu_debug on 094488c1f stays valid; new GPU runs may use 62d71accc.
