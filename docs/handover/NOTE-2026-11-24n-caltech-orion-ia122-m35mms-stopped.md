# NOTE 2026-11-24n Caltech -> orion: IA12.2 B m35mms stopped (handover to fresh rc3 m35mms)

- rc3 m35mms (02-11) passed 200 cycles: CW 0, FATAL 0, dt 12.008.
- The IA12.2 B m35mms stopped cleanly 00:30 (step 7333072.73, rc 0). Last cycle 42649, t ~3.52e5 (40.7% of tlim).
- rsts 00017-00020 are kept in prod_ia122_1123/m35mms/rst_keep. 00020 is the final one, at cycle 42649.
- hpc-sm-02-03 is free and reserved for RSG with 02-17 + 02-16. Tests run there until the rst arrives.
- RSG transfer: the USER will open the MPCDF gate on Caltech, then we scp rsg.01707.rst from viper ourselves.

-- Caltech
