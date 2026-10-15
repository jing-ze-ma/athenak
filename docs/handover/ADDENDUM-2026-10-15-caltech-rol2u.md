# ADDENDUM to TASK-2026-10-15-caltech-rol2u (Raven pre-test results; fills "RESULTS_PLACEHOLDER" in section 0). No step of the TASK changes.

Raven A100 gpudev, 1 node x 4 A100, 1 MeshBlock per GPU, same code (1bd0c2ce a100 md5 8a1f8afd; FIN4 5c504e49 = d505d21f numerics;
5b7e1bc5 a100 481754ef). Harness /raven/ptmp/jinma/caltech_rol2u_1015/jobs/{xrst,xrol}.sbatch, logs/ and runs/ there.
Seeds: E1 = a 5b7e1bc5-written rst (cycle 30, dt_retry on); He giant = a FIN4 5c504e49-written rst (cycle 24, no dt_retry key).

| check | E1 (5b7e1bc5 rst, +20 cyc) | He giant (FIN4 rst, +16 cyc) |
|---|---|---|
| GATE-BASE: production binary P vs 1bd0c2ce Q (default keys) | rst data BITWISE, user+hydro hst 20 rows BITWISE | rst BITWISE, hst 16 rows BITWISE |
| P vs QD = Q + time/dt_retry=true (absent on the HG rst: added on restart, "dt_retry on" printed) | BITWISE | BITWISE |
| GATE-RESTART ROL2U C vs L1 -> `-r` -> L2 | BITWISE (rst + last 10 hst rows) | BITWISE (rst + last 8 rows) |
| GATE-RETRY ROL2U C vs F (forced, fac 1.0) | BITWISE, 1 `DT_RETRY ... reason=forced` | BITWISE, 1 forced line |
| FATAL / NC / floor clips | 0 / 0 / 0 all arms | 0 / 0 / 0 all arms |
| Picard mean P -> ROL2U | 3.00 -> 3.02 (early E1 state, already 3) | 6.34 -> 3.00 |
| ROL2U tol_rel: predicted solves, mean tol_eff | 40, 2.8e-5 | 32, 1.0e-4 (= the implicit_tol_rel_max cap on every solve) |
| `<par_end>` offset seed -> P -> Q/QD -> ROL2U | 37 384 -> 37 338 -> 37 778 -> 37 813 (F 37 965) | 40 926 -> 40 898 -> 41 404 / 41 741 -> 41 776 (F 41 922) |

Conclusions for Caltech:
1. FIN4 (d505d21f) and fin4-dtretry (5b7e1bc5) rst ARE readable by 1bd0c2ce, and BASE from them is BITWISE to the production binary,
   also with `time/dt_retry=true` added -> GATE-BASE is expected to PASS on H200 too (a failure there = platform/build issue: report).
2. ROL2U restart and forced retry are BITWISE on both cases -> GATE-RESTART / GATE-RETRY expected PASS.
3. Header growth from a FIN4 HG rst: +0.85 KB (FIN5.1 dump +0.5, dt_retry keys +0.34, ROL2U keys +0.04 since GetOrAdd defaults are
   already dumped). Your 200-Rsun rst header 41.6 KB -> ~42.5 KB: under the gate's 44 000 limit; a forced-retry test rst adds 0.15 KB.
4. Watch item (He giant): tol_eff sits at the cap 1e-4 (vs implicit_tol 1e-8) on every predicted solve; this is the regime where the
   window comparison vs BASE (and the optional cfl/2 arm H) matters most -> please run H for the He giant if any budget allows.
