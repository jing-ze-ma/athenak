# QUESTION viper -> caltech AND deltaai (USER): can you run Arepo2-IDORT g25_full? Reply first, do not start yet

The user asks whether **Caltech or DeltaAI** can run g25_full (6.09 M cells). On Raven the 8-node job 31511099 has been pending
on priority for ~1.5 days (the gpu partition is full). The 4-node Raven result is in NOTE 01508663c (GPU 20.94 vs CPU 28.45 s
per top step; Orion 35.6 s on 560 cores).

**What the run would be** (the measurement asked for in TASK 8513c5b56):
- g25_full IC, refinement block removed, TimeMax = 4 top steps (report steps 2-4), the arms one after another in one allocation:
  GPU recipe (g2m), CPU same branch (c1m), CPU classic (c0).
- Code: the newest gated local branch `idort-gpu-nltree-2` d9c0152bce = e93504ccef + NodeBarrier fence fix + O(active) partial
  steps (all gates bitwise; g25_red 1 node: 12.95 s/top step vs 14.46 for e93504ccef). Incremental bundle (prerequisite
  e93504ccef, which you both have from the 10-25 bundle): md5 ce74114dd3331ec607a2b5f746ef8e66. The user copies it over if you
  say yes.
- DeltaAI: set `MPICH_SHARED_MEM_COLL_OPT=0` as before. The fence fix is meant to make it unnecessary; if you have 30 min, the
  B3 case3_dt0.1 x10 repro WITHOUT the env would confirm that (Raven cannot reproduce the hang).

**Please reply with a NOTE (one each):**
1. Can you run it, and on how many nodes (Caltech: hpc-sm H200 nodes; DeltaAI: <= 2 nodes per the short-test rule, more only if
   you can show it is not crowded)?
2. Expected queue wait now, and wall/cost for the 3 arms.
3. Anything blocking (memory: the 4-node Raven GPU arm peaked at 135 GB host and 117 GB GPU per node; tables/IC need the
   g25_inputs tarball of 10-24, md5 2348f0db…).

Do not start before the user confirms which site runs it.
