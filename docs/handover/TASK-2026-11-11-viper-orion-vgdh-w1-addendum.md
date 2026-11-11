# TASK addendum 2026-11-11 viper -> ORION: vgdh-w1 -- sp binary and layout (read BEFORE starting vgdh-w1 / vgdh-geo sp parts)

1. Every **sp / pulse_u arm** (vgdh-w1 step 2, vgdh-geo step 2, any SP_FULL ladder) uses **fork/vgdlay-1110 453f879bd**.
   Bug fixed there: on the sp wedge the vet_gd half-range data vgd_hr had zero ghosts (never exchanged), and the plm-on-h /
   hc-plm face stencils read them at every MeshBlock boundary and at the periodic x3 wrap. pulse_u full C_128 CAND: shell-mean
   E differed 9.5e-3 between 8x8 and 4x4 lateral blocks (fixed: 4e-5); a single block also moved 4.7e-3 (periodic wrap).
   All earlier pulse_u CAND numbers (viper and Raven SH 1.58/1.22, viper 0.82/0.76, the lf arms) are INVALID.
   Production sp keys (no implicit_hr_recon_q = h) are bitwise unaffected (pulse_u C_128 tab, old vs fixed).
2. **cs arms may stay on 92d4621fd**: the fix is sp-only; cs keys-absent bitwise vs 92d4621f (pulse_cs C_64, 10 cycles).
3. **Fix ONE lateral MeshBlock layout in physical angle for all sp levels** (e.g. lateral meshblock = 1/4 of the angular
   extent at every level, or one block per level); report the layout.
Build line unchanged (NOTE-2026-11-08-orion-viper-consfix-ledger.md). Details: viper /viper/ptmp2/jinma/rc18_1110/VGDH.md
"Layout bug". NOTE back as part of NOTE vgdh-w1.
