# TASK addendum 2 (2026-11-11) viper -> ORION: vgdh-w1 -- add the constant-w arm wc5

viper (pulse_cs C R3 32-64-128, CPU): CAND 0.91 / 1.11; CAND-w1 (w = 1) 2.84 / 2.68; **CAND-wc5 (w = 0.5 on every face)
2.58 / 2.58**; shell-mean g (implicit_thin_gface_shell, fork/vgdhr-1110 13cba7d2f) 0.91; formal-solution central branch
(implicit_hr_fsc) 0.34. So a constant w is 2nd order; the O(dx) is the radial VARIATION of w = S(g). Please add to vgdh-w1:
- **CAND-wc5**: `implicit_thin_h0 = -1.0e6`, `implicit_thin_h1 = 1.0e6` -> pulse_cs C R3 32..256 (+ T 8..64), on 92d4621fd
  (or 13cba7d2f: keys absent bitwise); pulse_u SP_FULL on fork/vgdlay-1110 453f879bd (or 13cba7d2f), one fixed layout.
- Cartesian beam suite with wc5 on top of the Cartesian cand (cb / ws66 / shd / hcyl, score.py), next to the w1 row.
No need to run implicit_thin_gface_shell / implicit_hr_fsc (not fixes). NOTE back inside NOTE vgdh-w1.
