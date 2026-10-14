# NOTE Caltech -> Raven (10-14 16:03 PDT): E1W links submitted (GO 95197aff, user OK); He giant C2/C3 dropped

- E1W: e1w_1014/prod/submit_links.sh -> agcE1w 4537666 (link 1, pending) + 4537667/68 (afterany); 1 H200, QoS normal, 7 d, tlim 1e9,
  binary 6d4215c0, IN remap/e1w.athinput. Side-by-side reports with e1s will follow once it runs (Picard, BE/stage fallbacks,
  DT_RETRY, s/cycle at production cadence, lateral cells per H_p).
- He giant LARGE: C2/C3 BG jobs 4522986/7 removed (one had just started; no results). FIN5.1 build + gate51 on FIN5.1 in progress
  (separate NOTEs); production waits for the gate result and the user's OK.
- e1s health flag: by t = 2.389e7 s run.log.4513330 has 361 DT_RETRY lines and 207 NON-CONVERGED (Picard 200-pass) events, s/cycle
  3.0 (was 1.8 at 2.1e7), dt 627 s (was ~730). No FATAL. Rate rose from ~1/1000 cycles before 2.0e7 to ~25/1000 at 2.2e7 and higher since.
  Full table at the 2.42e7 report.
