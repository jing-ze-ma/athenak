# NOTE Caltech -> Raven: beam baseline WITHDRAWN, acknowledged (2026-10-10 ~18:00 PDT)

Stopped on our side. The analysis job was cancelled and nothing is running for this task.
Raw runs are kept in case they save time. 89 of the 90 finished with rc.txt, at Caltech
/resnick/groups/carnegie_poc/jingze/arms_1010/beam_baseline/{x1,x2}/c1e{2,3,4}/<beam>/<arm>/ (run.log, rc.txt, time.txt, bin/),
binary athena_cpu_built_in_pgens_e50ec7d6 (f684f0c0).
- The x1 set (45 runs) finished in 2.7 h.
- The x2 set was held up by ONE run: x2 / c_light 1e4 (c dt/dx ~3000) / ba0 / hrx0. It took ~1350 s/cycle and was at
  cycle 30 of 200 after 11.2 h; we stopped it. All other x2 runs finished.
  That stall is the main finding: hrx0 at high c dt/dx on ba0 at x2 resolution is effectively unusable.
Tables were not produced. Ask if you want the metrics extracted from these dirs.
