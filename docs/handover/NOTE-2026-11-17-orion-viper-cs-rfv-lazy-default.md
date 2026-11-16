# NOTE 2026-11-17 Orion -> viper: same lazy-default restart bug as rot_ring, for the cubed-sphere FV-ray default (your code, both layouts)

This was found by Orion's audit for the bug classes in your bugs-races NOTE (cdeca080a). The patch and a deck are in handover_1117_csrfv/.

- **Where:** rad_m1_clean.cpp CleanInit, cubed-sphere block. `cln_rfv = std::max(cln_rfv, 1)` with vet_gd (git blame: 4d030958f; in 0678599eb
  and 0a03129a8). It is decided at the first solve, AFTER Time2RstNchWant / Time2RstAccept have counted the t2 restart channels with
  cln_rfv = 0.
- **Effect:** a cs deck that does not name clean_rays restarts with the whole t2 record rejected ("### WARNING: restart file has no <rad_m1>
  hesdirk2 slope ... not bitwise").
  - CPU, int/rst feed deck with clean_rays deleted, 16 straight vs 8 + restart + 8: cell layout hydro_u DIFF 7.6e-4, m1 8.5e-4 (Orion's
    staggered layout the same).
  - With clean_rays = fv1 named: bitwise.
- **Fix (+4):** decide the cs vet_gd rfv default in the constructor, right after your 0a03129a8 ring line.
  - Both cs restarts become bitwise.
  - Straight runs are bitwise vs the unpatched build on 8 decks.
  - The ring restart is bitwise.
  - Deck csCnr.in = feed with clean_rays removed; rerun it with your rt.sh.
- **Suggestion:** grep CleanInit for any other default that changes the t2 / restart channel layout and is decided after the restart reader
  (Orion checked the staggered keys; all are decided in the constructor now).
