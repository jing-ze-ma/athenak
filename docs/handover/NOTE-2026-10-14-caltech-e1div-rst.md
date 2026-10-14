# NOTE Caltech -> Raven DR: E1 restart arms from rst 00022 (TASK e1div-rst, part 2) -- QUEUED 10-13 19:20 PDT

- Dir e1div_rst/{R0,R1}; rst 00022 hard-linked (prod untouched). Binary dad00933 (ff7b7f11), launch e1_1013/e1launch.sh, 1 node x 2 H200.
- Keys exactly as the TASK (last_time 1.1e7 for outputs 1-5, output5/dt=1e30, output2/3 dt 2.5e4, -t 05:40:00);
  R0 tlim 1.26e7; R1 + time/cfl_number=0.15, tlim 1.30e7.
- Jobs: **R0 4483862, R1 4483863** (6 h each, PENDING Priority). Report + profiles follow when they end.
