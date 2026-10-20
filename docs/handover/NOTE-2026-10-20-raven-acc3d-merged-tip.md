# NOTE 2026-10-20 raven -> caltech: acc3d — switch to the merged production tip at the next link boundary (results identical)

USER rule (10-20): all production updates live on ONE merged branch, prodmerge-1019b (it becomes fin51-r-fix after its gates).
Raven checked: prodmerge-1019b @ 13206541 gives BITWISE the same accretor results as acctrc-1019c 5f29b077 on deck v11amwtr
(static C/C2/R, stream S/Sm/SR, forced retry SF: all hst incl. user/trc, bins and rst data incl. the history block identical).
Action: keep the running/queued chain (5086441-49) as is; build origin/prodmerge-1019b at SHA 13206541 (NOT the branch head, which
has an ungated cs-s7 merge on top) H200 nofma, report md5, and switch BIN at the next link boundary (no method change). Raven's
own fresh copy now uses 13206541 (A100 md5 29caa1ab).
