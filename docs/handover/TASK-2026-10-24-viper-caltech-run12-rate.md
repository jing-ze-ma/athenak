# TASK viper -> caltech (USER): report the run12 throughput (read-only, no change to the run)

The user wants to know which machine runs run12 faster. Please reply in a NOTE with:
1. **Static phase:** the wall time and final cycle (from `chain.txt` / the link log: the "static rc 0 wall ... cycle=" line).
2. **Stream phase now:** the last `elapsed= cycle= time= dt=` line of `run12/stream/run.log`, plus the cycle at which the stream started.
3. **Rates:** s per cycle in both phases, and simulated time per wall hour in the stream phase.
4. **Other runs:** anything else sharing the node or GPUs during these phases.

viper reference (link 1, job 12327914, 4 nodes x 2 MI300A = 8 GPUs, 8 ranks):
- **Static:** 4343 cycles in 3639 s, 0.84 s/cycle.
- **Stream:** elapsed 11284 s, cycle 17100, time 0.09809, dt 3.0e-6.
  - About 0.88 s/cycle if the cycle count continues from 4343.
  - About 0.0154 time units per wall hour.

Read-only: do not stop, restart or retune run12.
