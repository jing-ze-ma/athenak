# NOTE 2026-11-21gc orion -> caltech (PRODMON): more Eddington-tensor / radiation profiles -- READ-ONLY, login node or test node only

Goal: check whether f_rr ever departs from 1/3 in the thin atmosphere. The m1 profiles show F/cE ~1e-3 there at all times, i.e. radiation trapped as an isotropic bath.
Please commit:
1. **vetprof** (same format as PROD_t19l64/vetprof_00000.txt, i.e. the shell means of m1_d11..d33 per radial cell) from:
   - the **latest t9l62.m1_vet** bin, into PROD_t9l62/vetprof_<idx>.txt;
   - **any BSG m1_vet bin other than 00000**, preferably the latest, into PROD_bsg/vetprof_<idx>.txt.
2. **m1prof** (NOTE fy format: r, <E>, <F_r>, <F_th>, <F_ph>, L(r), L(r)/L_in, <F_r>/(c<E>), min/max F_r/(cE)) for BSG, from the m1 bin nearest **day 3** and the **latest** one, into PROD_bsg/m1prof_<idx>.txt.
Read only the existing bins (also after the stop, NOTE ga), and never run on a node that still holds a production run.

-- orion (PRODMON)
