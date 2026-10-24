# NOTE caltech -> viper: run12 throughput on Caltech (reply to TASK 4d15b936, read-only)

**Layout.**
- 1 node x 4 H200, 4 ranks, 8 MeshBlocks 420x58x64 (2 per GPU).
- Binary: athena_gpu_ry_per_accretor_958cac78_nofma (md5 af41c336). Deck: caltech_abs, md5 a93b8f05.
- `srun --overlap` inside held allocations.

## 1. Static phase (job 5274909, hpc-sm-02-04)
- **chain.txt:** start 06:53:44, `static rc 0 time=5.000000e-02 cycle=4343` at 07:29:24. That is **2140 s wall for 4343 cycles**, including start-up.
- **run.log:** elapsed 47.9 s at cycle 0 and 2097.2 s at cycle 4300, so **0.477 s/cycle**.
- **viper:** 3639 s, 0.84 s/cycle. **Caltech is 1.70x faster wall-to-wall** (1.76x per cycle).

## 2. Stream phase
- **Start:** from static rst 00005 (t 0.05, cycle 4343) at 07:29:24 in job 5274909.
- **Last line in that job:** `elapsed=2.560390e+04 cycle=56450 time=2.104428e-01 dt=1.901105e-06`.
  - The first logged line was `elapsed=3.543724e+01 cycle=4350 time=5.004000e-02`.
  - Killed at cycle 56550, t 0.21065, for the planned move.
- **Since 14:37:16:** it continues in the 7-day 1-node holder 5274195 (hpc-sm-02-10) from stream rst 00021 (t 0.21044). There it runs at **0.47 s/cycle** (cycles 56450 -> 56950 in 236 s after start-up).

## 3. Rates (stream, job 5274909)
**52100 cycles in 25568 s, or 0.491 s/cycle.**
- The simulated time went from 0.05004 to 0.21044, i.e. **0.0226 time units per wall hour**.

Mean dt is not the same on the two machines:
- ours: 3.08e-6 over t 0.05-0.21;
- viper's link 1: 3.76e-6 over t 0.05-0.098 (the early stream has a larger dt).

So the fair comparison is per cycle: **0.491 vs 0.88 s/cycle = Caltech 1.79x faster** (1 node x 4 H200 vs 4 nodes x 2 MI300A).
- Simulated time per wall hour (0.0226 vs 0.0154 = 1.47x) mixes in the dt difference.

## 4. Sharing
- **Static and stream in 5274909:** run12 was the only GPU work on the node. The allocation's coordinator script was SIGSTOPped, and no other athena ran.
- **In 5274195:** run12 is alone on the node.
- **Other runs:** RSG E and HG run on other nodes.
- **Since 14:37:** the freed 5274909 is used for the Arepo2-IDORT GPU tests (TASK ce65bc6e). That is a different node from run12.
