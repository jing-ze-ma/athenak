# NOTE 2026-10-22 viper -> caltech: merge coordination on fin51-r-fix (re 9c432f17)

Both sites are about to put a branch on top of fin51-r-fix c92d5b03:
- Caltech: faceopac-prod-1022 = c92d5b03 + de0ee02b (kirchhoff) + d33f5f59 (harmonic), opt-in.
- Viper (USER 10-22): prod-extrap2-1022 = c92d5b03 + <rad_m1>/time2_vet_col = extrap2 (opt-in; one closure build per step + f_K, q
  extrapolated from the last two solved states; cures the hesdirk2 cycle-to-cycle swing in stiff thin cells that `predict` produces on
  the accretor stream phase, +7.6 % cost vs rebuild +41 %; key-off bitwise vs c92d5b03). Being gated now.
Proposal to avoid two diverging fast-forwards: Caltech fast-forwards fin51-r-fix to faceopac-prod-1022 when its gates pass and posts the
hash; viper then rebases prod-extrap2-1022 onto that tip, re-gates (key-off bitwise vs your tip, he + accretor), and fast-forwards
fin51-r-fix again. Viper will not move fin51-r-fix before your face-opacity fast-forward NOTE.
Accretor/BSG harmonic A/B (your items 1-3): starting on viper.
