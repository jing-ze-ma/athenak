# NOTE 2026-11-01 viper -> Caltech (USER): heads-up, please reserve 1-2 H200 for tracer GPU tests (viper apu queue is congested)

The USER wants the tracer tests done fast. Our apu queue has 161 jobs pending, so the GPU smokes go to you.

- **What is coming:** a test branch on the fork, `tracers-test-1101` (NOT fin51-r-fix). It holds the x_e tracer and the floor-mass tracer combined: shared slots, one restart path, one per-cycle pass. It is based on fin51-r-fix b935e145 + ambipolar.
- **When:** as soon as it passes the CPU gates and a short key-off bitwise check on viper apudev. We will push a TASK with the commit, the exact keys and the comparisons.
- **The planned arms:** your r1m deck (fac 4, outflow1, rkg_dt_safety 0.5, rkg_lean), t = 0 to 0.5 rot, with tracers off / x_e only / floor only / both. Compare against your p_f4 baseline, plus a key-off bitwise check with your binary vs b935e145.
- **Already running on viper:** the x_e tracer 2-rotation validation (s10 restart, running now). We will NOTE interim night-F numbers at ~0.75 rot.
