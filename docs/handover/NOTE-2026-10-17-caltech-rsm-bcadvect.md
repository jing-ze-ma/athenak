# NOTE 2026-10-17 caltech -> raven: rs_m decks, implicit_bc_advect = false removed (answer to TASK 6ef8f0de)

Done, deck-only. Tree: /resnick/groups/carnegie_poc/jingze/wavefix_1013/caltech_wavefix/runs/cart/.
- The single line `implicit_bc_advect = false` was deleted from `<rad_m1>` in rs_m1.2/FIN3_n512, rs_m2/FIN3_n512,
  rs_m2/FIN3_n1024, rs_m3/FIN3_n512 and rs_m3/FIN3_n1024 (`in.athinput`). Each removed exactly 1 line, and the key no longer
  appears in any of them. Copies of the old decks: `in.athinput.bak_1017`.
- Generator templates: none in our tree carries this key. run_wavefix.sh and wv_one.sh only set the different key
  `implicit_bc_advect_wall=false` (bitZ arm), and those are left unchanged.
- Not rerun: no Caltech numbers for the five cases. Production is untouched.
