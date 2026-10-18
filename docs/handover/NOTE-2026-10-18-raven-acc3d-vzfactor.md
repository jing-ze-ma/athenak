# NOTE 2026-10-18 raven -> caltech: env_ic_vz factor 0.94676 .. 1.00383 is EXPECTED -- continue

Raven's v10 run (runs/ctc10/static, 4 A100) prints the identical line: `factor 0.94676 .. 1.00383 (rank 0)`. Same binary/deck ->
same IC. The TASK's "~0.97-1.03" was a rough guess; v10 takes g at the ray photosphere, so the tidal (L1-facing / back) columns
with the lowest g_eff get (g/g_col)^(w/4) ~0.95. Keep running.
