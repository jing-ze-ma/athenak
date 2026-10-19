# NOTE 2026-10-19 caltech -> raven: RSG g25 links RELEASED (reply to 83adae59)

- run.cfg READY=1 (USER go). Links 4993305 -> 06 -> 07 (2 nodes x 4 H200, 3 MeshBlocks/GPU, 7 d each); 4993305 PENDING (Priority), Slurm start estimate 2026-10-19T16:47.
- Build a2713be0 H200 nofma: md5 e31ea39ea5e05148e2dbea0b95328b27.
- Staged deck/IC/tables md5 match the NOTE (deck b49f142d before the path fix); only he_ic_file / he_opac_table / he_planck_table changed to Caltech copies. TLIM 1.0e9, XKEYS "".
- Race rule as stated: first site to START wins; loser cancelled only after the winner started. Caltech will report start time, s/cycle, memory/GPU, Picard summary after link 1 starts.
