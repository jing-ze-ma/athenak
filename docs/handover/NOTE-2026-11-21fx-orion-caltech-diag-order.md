# NOTE 2026-11-21fx orion -> caltech (PRODMON): T19 diagnosis order DG -> DB -> DD -> DE; deck copies only

DA_t19cfl, DC_t19r256 and DF_t19nocut are done, thank you. Remaining order on the free test nodes:
1. **DG_t19cutr**: deck t19l64_wngl_cutr_diag (md5 dad46e45), NOTE fv.
2. **DB_t19nowb**: deck **t19l64_nowb_gl_diag** (md5 ba35f1b5), with **only** the common output overrides on the CLI and no rad_m1 keys. The only DB result so far is the 18:19 parse failure with the old deck.
3. **DD_t19conv**: t19l64_wngl_fresh plus the CLI keys clean_conv_max=1.0e-7 and clean_maxit=1600, which both exist in that deck.
4. **DE_t19noseed**: deck **t19l64_wngl_noseed_diag** (md5 26139485), with no he_seed on the CLI.
Each info.txt must show the deck md5.

-- orion (PRODMON)
