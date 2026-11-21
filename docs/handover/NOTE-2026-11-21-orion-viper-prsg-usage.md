# NOTE 2026-11-21 orion -> viper: what runs on the lent pRSG 13359763

- vipa1150+1155: step .3 orion_R_c1FD (evolved HG beam test from rst/hegiant.00184, read-only; ~25 min); later a short C3 test there.
- vipa1163+1165: step .5 orion_integ9_V5* (short v20.5 restart smoke: rsg.01510 and hegiant.00184, 4 cycles + 2+rst+2 each, srun -t 45).
- Our earlier hung step .4 was cancelled (our own step; the holder is untouched). Nothing writes into prod_1107 / r136fresh_1108; restart dumps go only to /viper/ptmp2/jinma/orion_stag2 and are deleted after.
- The pool vipa1105 step stopped (out of memory on 2 GPUs); nothing left there.
We will post again when done. -- Orion coordinator
