# NOTE 2026-11-23zz17 Caltech -> orion: m35mms production switched to IA12.2 B (USER OK on Caltech)

- IA11 m35mms stopped cleanly at 22:04 (cycle 36376). Final rst m35mms.00015.rst is hard-linked in prod_ia11_1123/m35mms/rst_keep and in the new dir.
- IA12.2 B continuation started ~22:05 on hpc-sm-02-03 (step 7333072.73), dir prod_ia122_1123/m35mms/.
  - Binary athena_ia122_465aff485 (md5 cfa3dc77…).
  - Restart -r 00015 with -i docs/handover/ia122_decks/notile/IA122_m35mms.athinput (names every IA12.2 key) + rad_m1/vet_gd_wall_sym=true rad_m1/vet_gd_theta_seam_sym=true. Both keys are echoed "added from the command line (= true)".
  - The deck keeps hydro_w/m1 bins every 1800 s, rst every 21600 s, KEEPRST 3. Launcher: prod_ia122_1123/launch_ia122.sh.
- First cycles: 0 FATAL, 0 CW. dt 12.006 (IA11's last was ~6.6).
  - ~1.35 s/cyc in the first window, including the startup window, against IA11's 1.26.
  - So ~1.7x more model time per wall second. A clean s/cyc follows in the PROD table.
- The IA11 PROD-table poster still reads prod_ia11_1123; we'll point it at the new dir.

-- Caltech
