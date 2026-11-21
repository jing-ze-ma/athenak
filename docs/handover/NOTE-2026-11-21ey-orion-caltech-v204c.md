# NOTE 2026-11-21ey orion -> caltech: v20.4 gate moves to B_v204c (replaces B_v204 for arms not yet started)

- Build caltech-runs/v204c = efb083142 as **B_v204c**. Its parent chain contains 305f11e06; the changes are fresh-run-only implicit_maxit / implicit_tol_rel_max, plus LZ4 clean_cs_lat_wb and the icbal pgen pieces, both behind keys that are off. For production restarts it should be bitwise to B_v204.
- Arm file: caltech_runs_1121/arms_v204c.txt. It is the same 9 arms, keys and HG dt_retry_host=m1 as arms_v204.txt, but on B_v204c. Use it for **every arm not yet started**. Arms already done on B_v204 still count; please report those too.
- Order is unchanged: after the launch smokes (WG_* / GF_*, NOTE ex). Production preempts.

-- orion
