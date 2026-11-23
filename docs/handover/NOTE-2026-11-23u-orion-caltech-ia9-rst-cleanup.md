# NOTE 2026-11-23u orion -> Caltech: USER OK to delete the IA9 t=0 rst files

USER (11-23 ~13:25 PST): ok.
Please delete only the t=0 rst of each IA9 run, which can be regenerated from the deck:
- prod_ia9_1123/bsg/rst/bsg3d.00000.rst (21.8 GB)
- prod_ia9_1123/m35mms/rst/m35mms.00000.rst (28.5 GB)
- prod_ia9_1123/t9l62/rst/t9l62.00000.rst (86.3 GB)

KEEP the final rst of each (bsg3d.00001, m35mms.00001, t9l62.00001). t9l62.00001 is the start state for the diagnostic restart in 23q (arm A).
Please also keep the IA9 hst and logs. Nothing else changes.

-- orion (coordinator)
