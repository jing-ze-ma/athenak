# NOTE 2026-10-09 Caltech -> viper: ALL long runs on HOLD until the updated RT scheme (user 10-09)

User: "I think we can wait before the updated rt scheme. hold all long runs including bsg".
- Seeded BSG true repro (prod_seed): chain 4294781-83 cancelled at 16:10 PDT, 2.5 min into link 1 (cycle 50, t 4.4e3 s);
  only rst bsg3d.00000.rst (t = 0). CANCEL file in bsg_1009/prod_seed. Input (he_seed 1e-6 block, md5 a3d332ab) and smoke
  results kept; it will start fresh on the updated scheme when you send the binary/keys.
- Unseeded BSG (prod_noseed, 6.45 d, symmetric reference) and its page https://claude.ai/artifact/BQ8LZRU8uV1Zw83GtZQQxK stay.
- He giant fresh N897: already stopped (rst 1.0 d, NOTE f4b5c36f). Old remapped N897 continuation: stopped at 23.50 d (rst 00048).
- Accretor RHD (accretor-rhd-1008, merged mem-1009): option (a) transverse stream T profile implemented and gated, but the
  injected-radiation stream still diverges where the stream head meets the atmosphere (implicit M1 + vet_col in a radial
  column with outward free-streaming photosphere and an inward-shining opaque stream head; design doc sect. 11.3, dumps
  rhd/s5/crash_prof). Waiting for the updated scheme as well; no accretor jobs running.
- Short jobs continue: the xthinfix battery (TASK-2026-10-09-caltech-xthinfix-tests) is in progress; profiling done (7a20aeec).
  Caltech can take further SHORT tests; long productions resume when you send the updated-scheme task.
