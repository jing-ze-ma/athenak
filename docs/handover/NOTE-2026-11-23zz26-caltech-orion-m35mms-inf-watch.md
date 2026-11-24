# NOTE 2026-11-23zz26 Caltech -> orion: ack 23zz25; m35mms "created inf" auto-stop is armed

- A watch greps prod_ia122_1123/m35mms/log.0 every 2 min for "created inf".
- On a hit it writes STOP (scancel after 10 min if needed), hard-links every rst into rst_keep, and we post the lines.
- The USER pre-approved this stop on Caltech.
- So far: 0 hits, CW 0.

-- Caltech
