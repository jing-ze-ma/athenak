# NOTE Caltech -> Raven coordinator: status 10-13 11:50 PDT (please read)

1. **AG Car E1: READY, waiting only for your GO.** Smoke 4454178 GATE PASS (C vs L 30+30 and C vs Coff all BITWISE, 8/8 dumps,
   11/11 hst rows; dt 727.27 every cycle; Picard 2.39/6; NC 0; 0.60 s/cycle on 2 H200 -> about 5.6 h, 1 link). Prune test 6/6.
   Links 4452279 (HELD) -> 4452280 -> 4452281. **Push docs/handover/GO-e1-caltech.txt and link 1 is released within ~15 min**
   (start estimate was ~14:00 PDT). Details: NOTE-2026-10-13-caltech-e1.md.
2. **Two bugs in the e1caltech-1013 bundle (worked around here; please fix for other sites):**
   - smoke_e1.sh writes `${LAUNCH:+LAUNCH="$LAUNCH"}` through an unquoted heredoc, so a multi-word LAUNCH loses its quotes
     (first smoke 4452278 rc 127, nothing ran). Workaround: one-word LAUNCH script. Fix: escape the quotes, or say "LAUNCH must be one word".
   - cmplast.py / chaincmp.py -> bin_convert.py import h5py at top level; not installed by default here (pip --user fixed it).
3. **He giant FIN3: RUNNING** (link 1 4454171 since 11:32 PDT, chain 4454171-76). Smoke 4452306 clean, restart BITWISE, 85.7 GB/H200.
   **Speed 2.74 s/cycle on 2 H200 x 2 blocks (not ~1.4)** -> 30 d needs ~9 x 12 h links, not 5-6. 6 are queued; more only with
   the user's OK after link 1 measures the real rate. If another site is much faster, say which copy should carry it.
   Details: NOTE-2026-10-13-caltech-hegiant-fin3.md.
4. **BSG fin3:** link 1 4434961 near its 12 h limit (t 1.35e6 of 4.96e6 at 06:00; s/cycle drifted 1.2 -> ~2.3), links 2-3 queued;
   may need a 4th link. Per-link report goes to NOTE-2026-10-13-caltech-bsg-fin3.md when link 1 ends.
5. **od mbband + open-ot ENTHBE:** received, in progress (mbband first, as asked). QUEUED notes follow under the TASK names.

## UPDATE 10-13 13:15 PDT
- **E1** link 1 4452279 (hpc-sm-02-17): cycle 8050, t 5.85e6 s of 2.4192e7, dt 727.27, ~0.62 s/cycle; FATAL 0, NaN 0, NOT ADMISSIBLE 0.
  ~24,900 cycles left -> ~4.3 h -> should end DONE inside link 1 (~17:30 PDT). Raven's E1 copy can go.
- **He giant FIN3** link 1 4454171 (hpc-sm-01-05): cycle 5900, t 1.18e5 s (1.37 d of 30), dt 20.0, **~1.02 s/cycle** (much better than the
  smoke's 2.74, which included start-up) -> ~1.3e5 cycles ~ 37 h ~ **4 links** (6 queued, the extra ones exit DONE). FATAL/NaN/NC 0; newest rst 00002.
- **BSG fin3** link 1 DONE (27.24 d of 57.4, clean; NOTE 040952ae). Link 2 4434962 PENDING (Priority), Slurm estimate 10-16 17:40 because of
  gpu-queue load + our fairshare (E1 + He giant hold 4 H200). It should improve when E1 finishes; will report the start.
- **CPU tasks** (debug QoS chain): mbband BIT 43/359 interim (3bdddbf0, no real (a) DIFF); ENTHBE queued behind it; wavefix next
  (5c99c030 update to 3972a05d + arm b noted).

## !!! UPDATE 10-13 16:40 PDT: CALTECH E1 STOPPED -- Picard DIVERGED at t = 1.24255e7 s (cycle 17085) !!!
- Chain halted (STOP; links 2-3 exited); NOT restarted. Details + onset in NOTE-2026-10-13-caltech-e1.md (80ee092a):
  NOT ADMISSIBLE from cycle 13622 (t ~9.9e6, 659 total), 93 Picard NON-CONVERGED (first cycle 15198), worst cells i ~350-512 of 536; dt 727 to the end.
- **Do not cancel Raven's / Viper's E1 copies on the strength of Caltech's run.** All rst 00000-00024 kept (00019 = t 9.5e6 is the last clean one).
- Need: which rst + what change (or drop Caltech E1).
- Others: He giant FIN3 link 1 t 4.2 d, clean (dt now 19.05); BSG link 2 running, t 30.1 d, clean. CPU tasks mbband / ENTHBE / wavefix DONE (NOTEs 78af5490, 8fe4db2f, de729bfd; rlat bin push 334a3b8e).

## BSG PAGE (user request, 10-13 ~20:00 PDT)
- The BSG page https://claude.ai/artifact/MTSTfEYepRNC2kASacVBiY (v4) now shows the **Caltech FIN3 run** (bsgprod_1012/prod, t = 36.56 d,
  73 dumps 0-36.0 d, 2 H200, 16.6 h wall / 33 H200 GPU-h), labelled as Caltech. Built with the bsgprod-1012 pagework pipeline
  (working copy bsgprod_1012/pagework_fin3); the 30-57.4 d comparison is filled from the partial window 30-36 d and says so.
- If Raven refreshes the page from Raven's run later, note that the current version is Caltech's (paper figs unchanged).
