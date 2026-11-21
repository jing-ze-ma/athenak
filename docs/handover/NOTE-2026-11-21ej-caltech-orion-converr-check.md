# NOTE 2026-11-21ej Caltech -> orion: CONVERR check (answer to 604b04b06)

1. **Yes, the CET1_rsg rerun finished.** It ran 15:22-15:23 PST 11-20 on hpc-sm-02-04 (holder 6133592), binary athena_v202_f4cefb181
   md5 54c3df35aeb65bced81ad681e47a8020: rc=0, wall 83 s, 1 cycle, bins hydro_w/m1 06326 + m1_vet 01266.
   Our 16:00 line calling it "held" was stale. Its info.txt, res.txt and a 200-line run.log tail are now committed to caltech_results/CET1_rsg/.
   The earlier OOM attempt's tail was renamed OOM_attempt_1503_tail100.txt.
2. **Yes, CEP0_rsg is the 15:30 result (d73022852).** It ran 15:24-15:25 on 02-04 with the same binary: rc=0, wall 83 s, bins present.
So the CONVERR remainder is CEP0_hg + the CEPc/CET3c x4, as you say.

-- Caltech
