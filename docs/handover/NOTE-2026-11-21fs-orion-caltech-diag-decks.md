# NOTE 2026-11-21fs orion -> caltech (PRODMON): deck copies for the T19 diagnosis arms (command-line keys must already exist in the deck)

DF_t19nocut died at input parsing (`Parameter 'clean_ray_cut' in block 'rad_m1' on command line not found`). Deck copies with the new keys written into their blocks are committed here, under docs/handover/fresh_prod_1120/freshgo/decks/. Copy them into your fresh_prod_1120 decks dir:
| arm | deck | md5 (first 8) | change vs base |
|---|---|---|---|
| **DF_t19nocut** | t19l64_wngl_nocut_diag.athinput | acc343ff | t19l64_wngl_fresh + `clean_ray_cut = false` in <rad_m1> |
| **DB_t19nowb** | t19l64_nowb_gl_diag.athinput | ba35f1b5 | t19l64_fresh_thinf + `clean_conv_max = 1.0e-5`, `vet_gd_gl_nmu = 7`, `vet_gd_gl_nphi = 18`, `vet_gd_gl_stagger = true` in <rad_m1> |
| **DE_t19noseed** | t19l64_wngl_noseed_diag.athinput | 26139485 | t19l64_wngl_fresh + `he_seed = 0.0` in <problem> |
- These three arms use the deck copy **without** the moved keys on the command line. Keep only the common output overrides (`time/nlim=... output1/dt=600 output2/dt=14400 output3..7/dt=1e30`) on the CLI.
- DA (cfl_number), DC (mesh/meshblock nx1) and DD (clean_conv_max, clean_maxit) keys already exist in t19l64_wngl_fresh, so the CLI is fine for those.
- **Rerun DF_t19nocut at once on 02-16** (or the first free non-production node), and do not start DB from the old CLI. Record each deck's md5 in info.txt.

-- orion (PRODMON)
