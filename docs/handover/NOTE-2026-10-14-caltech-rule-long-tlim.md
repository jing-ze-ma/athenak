# NOTE Caltech -> Raven: USER RULE (10-14): make tlim very long, and walltime as long as allowed

USER (10-14): "yes we should make tlim as long as possible", "also raise tlim for e1", "let raven know that rule to make tlim very long".

Rule for production runs on every machine:
- **tlim very long** (e.g. 1e9 s), so a run is limited by its walltime links, not by a planned end time. Exception: a run with a fixed comparison span. The Linhao BSG keeps 4.96e6 (57.4 d, paper comparison); the user chose "He giant only" for that chain.
- **Walltime as long as the QoS allows.** On Caltech that is 7 d (normal). QoS long allows 14 d but only 1 job per user at lowest priority. athena -t = walltime minus a margin, so the final rst is written.
- **How to raise tlim:** set `time/tlim` in the deck (or in a read-only production deck copy), or on the command line on restart, since the key always exists. Also set the link script's TLIM (DONE check) to the same value. A DONE file left by the old tlim must be removed before the next link.

Applied on Caltech: He giant FIN4 chain (TLIM 1e9; 7-d links 4513229, 4513230 after 4489939), E1 seed 1e-3 (e1s-1014: production deck tlim 1e9, longest walltime). Details are in NOTE-2026-10-14-caltech-wide4-trials.md.
Please apply this to Raven/Viper/DeltaAI productions and to future TASK docs.
