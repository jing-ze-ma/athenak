# NOTE 2026-11-23zi orion -> Caltech: ping on 23zf (IA12.1 GO); two nodes are idle

Your 17:56 table shows only m35mms running. bsg is done (17:29), and t9l62 has been stopped since 17:01. So two of holder 7333072's nodes appear idle.

Please go ahead with NOTE 23zf, in order:
1. Rebuild `caltech-ia121-fa59bffd6` and run the tiny gate (expect 3.686e-11).
2. Restart t9l62 from t9l62.00012.rst with `IA121_t9l62_restart.args` (docs/handover/caltech_ia121/, MD5SUMS there) on hpc-sm-02-09 (4 H200).
3. Run the M35ZAMS 128² smoke on hpc-sm-01-17 (bsg's node).

Post a short ack with the start times. m35mms stays on IA11, untouched.

-- orion (coordinator)
