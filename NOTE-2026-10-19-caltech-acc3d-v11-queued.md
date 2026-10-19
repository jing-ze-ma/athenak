# NOTE 2026-10-19 caltech -> raven: 3-D accretor v11 queued on H200 (reply to TASK 98f0d042)

- Dir /resnick/groups/carnegie_poc/jingze/acc3d_v11_1019; deck = acc3d_prod_v11.athinput.in with @FILES@ -> files/ (column/tables md5 5b813886 / d81b7a46 / 99d668f6, unchanged from v3). Layout as v3: 2 H200, 8 MeshBlocks.
- Build 5007428: acc3d-held-1018 @ 23e06019, ry_per_accretor, nofma (md5 to follow).
- Check 5007429 (afterok build): static phase to t 0.002 C, C2, R from the t 0.001 rst; hst every cycle + final bins bitwise; exits 1 on any DIFF.
- Links 5007430 -> 31 -> 32 -> 33 -> 34 -> 35 (afterok chain from the check; 12 h each; acc3d_chain.sh unchanged: static to 0.05, then stream on), FRESH from t=0.
- Both sites kept (USER 10-18). RSG g25 jobs unchanged.
